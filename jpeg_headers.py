from collections import defaultdict
from typing import List, Tuple, Dict, BinaryIO


# Zig-zag scan order (# https://www.w3.org/Graphics/JPEG/itu-t81.pdf, Figure A.6)
zigzag = [(0, 0), (0, 1), (1, 0), (2, 0), (1, 1), (0, 2), (0, 3), (1, 2),
          (2, 1), (3, 0), (4, 0), (3, 1), (2, 2), (1, 3), (0, 4), (0, 5),
          (1, 4), (2, 3), (3, 2), (4, 1), (5, 0), (6, 0), (5, 1), (4, 2),
          (3, 3), (2, 4), (1, 5), (0, 6), (0, 7), (1, 6), (2, 5), (3, 4),
          (4, 3), (5, 2), (6, 1), (7, 0), (7, 1), (6, 2), (5, 3), (4, 4),
          (3, 5), (2, 6), (1, 7), (2, 7), (3, 6), (4, 5), (5, 4), (6, 3),
          (7, 2), (7, 3), (6, 4), (5, 5), (4, 6), (3, 7), (4, 7), (5, 6),
          (6, 5), (7, 4), (7, 5), (6, 6), (5, 7), (6, 7), (7, 6), (7, 7)]


class HuffmanTable:
    table_class: int
    dest_id: int
    counts: List[int]
    huff_data: Dict[Tuple[int, int], int]


class QuantizationTable:
    precision: int
    dest_id: int
    table: List[List[int]]


class FrameComponent:
    identifier: int
    sampling_factor: int
    h_sampling_factor: int
    v_sampling_factor: int
    quant_table_dest: int


class StartOfFrame:
    precision: int
    num_lines: int
    samples_per_line: int
    components: List[FrameComponent]


class ScanComponent:
    selector: int
    dc_table: int
    ac_table: int


class StartOfScan:
    components: List[ScanComponent]
    spectral_selection_range: Tuple[int, int]
    successive_approximation: int


def define_huffman_table(file: BinaryIO) -> List[HuffmanTable]:
    # https://www.w3.org/Graphics/JPEG/itu-t81.pdf, B.2.4.2
    
    # Lh: Huffman table definition length – Specifies the length of all Huffman table parameters
    size = int.from_bytes(file.read(2), byteorder="big")

    huff_tables = []

    bytes_left = size - 2
    while bytes_left > 0:
        huff_table = HuffmanTable()
        huff_table.huff_data = {}

        table = int.from_bytes(file.read(1), byteorder="big")
        # Tc: Table class – 0 = DC table or lossless table, 1 = AC table
        huff_table.table_class = table >> 4

        # Th: Huffman table destination identifier – Specifies one of four possible destinations at the decoder into
        # which the Huffman table shall be installed
        huff_table.dest_id = table & 0b1111

        # Li: Number of Huffman codes of length i – Specifies the number of Huffman codes for each of the 16 possible
        # lengths allowed by this Specification. Li’s are the elements of the list BITS.
        huff_table.counts = [int.from_bytes(file.read(1), byteorder="big") for _ in range(16)]

        length_codes_map = defaultdict(list)
        # Incrementing code used to build Huffman map
        code = 0
        for i in range(16):
            # Go through number of codes in each length
            for j in range(huff_table.counts[i]):
                huff_byte = int.from_bytes(file.read(1), byteorder="big")
                huff_table.huff_data[(code, i+1)] = huff_byte
                length_codes_map[i+1].append(huff_byte)

                # increment code
                code += 1

            # Shift bits to left (increase length by 1)
            code <<= 1

        bytes_left -= 17 + sum(huff_table.counts)

        print(f"Huffman table length: {size}")
        print(f"Destination ID: {huff_table.dest_id}")
        print(f"Class = {huff_table.table_class} (" +
              ("DC / Lossless table" if huff_table.table_class == 0 else "AC table") + ")")

        for i in range(16):
            print(f"    Codes of length {i+1} bits ({huff_table.counts[i]} total): ", end="")
            for huff_byte in length_codes_map[i+1]:
                print(f"{huff_byte:02X} ", end="")
            print()

        print(f"Total number of codes: {sum(huff_table.counts)}")
        print()

        huff_tables.append(huff_table)

    return huff_tables


def define_quantization_table(file: BinaryIO) -> List[QuantizationTable]:
    # https://www.w3.org/Graphics/JPEG/itu-t81.pdf, B.2.4.1

    # Lq: Quantization table definition length – Specifies the length of all quantization table parameters
    size = int.from_bytes(file.read(2), byteorder="big")

    quant_tables = []

    bytes_left = size - 2
    while bytes_left > 0:
        quant_table = QuantizationTable()
        quant_table.table = {}

        # Create empty 8x8 table
        quant_table.table = [[0 for _ in range(8)] for _ in range(8)]

        temp = int.from_bytes(file.read(1), byteorder="big")
        # Pq: Quantization table element precision – Specifies the precision of the Qk values. Value 0 indicates 8-bit
        # Qk values; value 1 indicates 16-bit Qk values. Pq shall be zero for 8 bit sample precision P.
        quant_table.precision = temp >> 4
        element_bytes = 1 if quant_table.precision == 0 else 2

        # Tq: Quantization table destination identifier – Specifies one of four possible destinations at the decoder
        # into
        # which the quantization table shall be installed
        quant_table.dest_id = temp & 0b1111

        bytes_left -= 65 + (64 * quant_table.precision)

        for i in range(64):
            # Qk: Quantization table element – Specifies the kth element out of 64 elements, where k is the index in
            # the zigzag ordering of the DCT coefficients.
            # The quantization elements shall be specified in zig-zag scan order.
            element = int.from_bytes(file.read(element_bytes), byteorder="big")

            row, col = zigzag[i]
            quant_table.table[row][col] = element

        print(f"Quantization table length: {size}")
        print("Precision: " + ("8 bits" if quant_table.precision == 0 else "16 bits"))
        print("Destination ID: " + str(quant_table.dest_id) +
              (" (Luminance)" if quant_table.dest_id == 0 else " (Chrominance)"))
        for i in range(len(quant_table.table)):
            print(f"    DQT, Row #{i}: " + "".join(str(element).rjust(4) for element in quant_table.table[i]))
        print()

        quant_tables.append(quant_table)
    return quant_tables


def read_start_of_frame(file: BinaryIO) -> StartOfFrame:
    # https://www.w3.org/Graphics/JPEG/itu-t81.pdf, B.2.2
    sof = StartOfFrame()
    sof.components = []

    # Lf: Frame header length – Specifies the length of the frame header
    size = int.from_bytes(file.read(2), byteorder="big")

    # P: Sample precision – Specifies the precision in bits for the samples of the components in the frame.
    sof.precision = int.from_bytes(file.read(1), byteorder="big")

    # Y: Number of lines – Specifies the maximum number of lines in the source image. This shall be equal to the
    # number of lines in the component with the maximum number of vertical samples (see A.1.1). Value 0 indicates
    # that the number of lines shall be defined by the DNL marker and parameters at the end of the first scan (see
    # B.2.5).
    sof.num_lines = int.from_bytes(file.read(2), byteorder="big")

    # X: Number of samples per line – Specifies the maximum number of samples per line in the source image. This
    # shall be equal to the number of samples per line in the component with the maximum number of horizontal
    # samples (see A.1.1).
    sof.samples_per_line = int.from_bytes(file.read(2), byteorder="big")

    # Nf: Number of image components in frame – Specifies the number of source image components in the frame.
    # The value of Nf shall be equal to the number of sets of frame component specification parameters (Ci, Hi, Vi,
    # and Tqi) present in the frame header.
    sof.num_frame_components = int.from_bytes(file.read(1), byteorder="big")

    for i in range(sof.num_frame_components):
        component = FrameComponent()

        # Ci: Component identifier – Assigns a unique label to the ith component in the sequence of frame component
        # specification parameters. These values shall be used in the scan headers to identify the components in the
        # scan. The value of Ci shall be different from the values of C1 through Ci − 1.
        component.identifier = int.from_bytes(file.read(1), byteorder="big")

        component.sampling_factor = int.from_bytes(file.read(1), byteorder="big")

        # Hi: Horizontal sampling factor – Specifies the relationship between the component horizontal dimension
        # and maximum image dimension X (see A.1.1); also specifies the number of horizontal data units of component
        # Ci in each MCU, when more than one component is encoded in a scan.
        component.h_sampling_factor = component.sampling_factor >> 4

        # Vi: Vertical sampling factor – Specifies the relationship between the component vertical dimension and
        # maximum image dimension Y (see A.1.1); also specifies the number of vertical data units of component Ci in
        # each MCU, when more than one component is encoded in a scan.
        component.v_sampling_factor = component.sampling_factor & 0b1111

        # Tqi: Quantization table destination selector – Specifies one of four possible quantization table destinations
        # from which the quantization table to use for dequantization of DCT coefficients of component Ci is retrieved.
        # If the decoding process uses the dequantization procedure, this table shall have been installed in this
        # destination by the time the decoder is ready to decode the scan(s) containing component Ci. The destination
        # shall not be respecified, or its contents changed, until all scans containing Ci have been completed.
        component.quant_table_dest = int.from_bytes(file.read(1), byteorder="big")

        sof.components.append(component)

    print(f"Frame header length: {size}")
    print(f"Precision: {sof.precision}")
    print(f"Number of lines: {sof.num_lines}")
    print(f"Samples per line: {sof.samples_per_line}")
    print(f"Image size: {sof.samples_per_line} x {sof.num_lines}")

    print(f"Number of image components: {sof.num_frame_components}")
    for i, component in enumerate(sof.components):
        print(f"    Component {i+1}: ID=0x{component.identifier:X}, "
              f"Sampling factor=0x{component.sampling_factor:X}, "
              f"Vertical sampling factor=0x{component.v_sampling_factor:X}, "
              f"Horizontal sampling factor=0x{component.h_sampling_factor:X}, "
              f"Quantization table destination=0x{component.quant_table_dest}")
    print()

    return sof


def read_start_of_scan(file: BinaryIO) -> StartOfScan:
    # https://www.w3.org/Graphics/JPEG/itu-t81.pdf, B.2.3
    sos = StartOfScan()
    sos.components = []

    # Ls: Scan header length – Specifies the length of the scan header.
    size = int.from_bytes(file.read(2), byteorder="big")

    # Ns: Number of image components in scan – Specifies the number of source image components in the scan. The
    # value of Ns shall be equal to the number of sets of scan component specification parameters (Csj, Tdj, and Taj)
    # present in the scan header.
    sos.num_scan_components = int.from_bytes(file.read(1), byteorder="big")

    for i in range(sos.num_scan_components):
        component = ScanComponent()

        # Csj: Scan component selector – Selects which of the Nf image components specified in the frame parameters
        # shall be the jth component in the scan. Each Csj shall match one of the Ci values specified in the frame
        # header, and the ordering in the scan header shall follow the ordering in the frame header.
        component.selector = int.from_bytes(file.read(1), byteorder="big")

        temp = int.from_bytes(file.read(1), byteorder="big")

        # Tdj: DC entropy coding table destination selector – Specifies one of four possible DC entropy coding table
        # destinations from which the entropy table needed for decoding of the DC coefficients of component Csj is
        # retrieved. The DC entropy table shall have been installed in this destination (see B.2.4.2 and B.2.4.3) by the
        # time the decoder is ready to decode the current scan. This parameter specifies the entropy coding table
        # destination for the lossless processes.
        component.dc_table = temp >> 4

        # Taj: AC entropy coding table destination selector – Specifies one of four possible AC entropy coding table
        # destinations from which the entropy table needed for decoding of the AC coefficients of component Csj is
        # retrieved. The AC entropy table selected shall have been installed in this destination (see B.2.4.2 and
        # B.2.4.3) by the time the decoder is ready to decode the current scan. This parameter is zero for the
        # lossless processes.
        component.ac_table = temp & 0b1111

        sos.components.append(component)

    # Ss: Start of spectral or predictor selection – In the DCT modes of operation, this parameter specifies the first
    # DCT coefficient in each block in zig-zag order which shall be coded in the scan. This parameter shall be set to
    # zero for the sequential DCT processes. In the lossless mode of operations this parameter is used to select the
    # predictor
    spectral_selection_start = int.from_bytes(file.read(1), byteorder="big")

    # Se: End of spectral selection – Specifies the last DCT coefficient in each block in zig-zag order which shall be
    # coded in the scan. This parameter shall be set to 63 for the sequential DCT processes. In the lossless mode of
    # operations this parameter has no meaning. It shall be set to zero.
    spectral_selection_end = int.from_bytes(file.read(1), byteorder="big")

    sos.spectral_selection_range = (spectral_selection_start, spectral_selection_end)

    sos.successive_approximation = int.from_bytes(file.read(1), byteorder="big")

    print(f"Scan header length: {size}")
    print(f"Number of image components: {sos.num_scan_components}")
    for i, component in enumerate(sos.components):
        print(f"    Component {i+1}: selector=0x{component.selector:02X}, "
              f"table={component.dc_table}(DC),{component.ac_table}(AC)")
    print(f"Spectral selection: {sos.spectral_selection_range[0]} .. {sos.spectral_selection_range[1]}")
    print(f"Successive approximation: 0x{sos.successive_approximation:02X}")
    print()

    return sos
