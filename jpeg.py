import jpeg_headers
from typing import Dict, Tuple, List
from time import time
import math
import tkinter

# References
#   https://www.w3.org/Graphics/JPEG/itu-t81.pdf
#   https://impulseadventure.com

IMAGE_TO_OPEN = "test_images/lena.jpg"

start_time = time()


class PhotoDisplay:
    def __init__(self, tk, height, width, rgb_array: List[List[Tuple[int, int, int]]]):
        self.width = width
        self.height = height

        self.i = tkinter.PhotoImage(width=width, height=height)
        pixels = " ".join(("{" + " ".join(("#%02x%02x%02x" % rgb_array[color_row][color_col]
                                           for color_col in range(width))) + "}" for color_row in range(height)))
        self.i.put(pixels, (0, 0, self.width, self.height))

        self.canvas = tkinter.Canvas(tk, width=self.width, height=self.height)
        self.canvas.pack()
        self.canvas.create_image(0, 0, image=self.i, anchor=tkinter.NW)

        self.canvas.bind("<ButtonPress-1>", self.scroll_start)
        self.canvas.bind("<B1-Motion>", self.scroll_move)

    def scroll_start(self, event):
        self.canvas.scan_mark(event.x, event.y)

    def scroll_move(self, event):
        self.canvas.scan_dragto(event.x, event.y, gain=1)


def bit_from_bytearray(arr: bytearray, bit_idx: int, order: str) -> int:
    if order == "little":
        return (arr[bit_idx // 8] & (0b1 << (bit_idx % 8))) >> (bit_idx % 8)
    else:
        return (arr[bit_idx // 8] & (0b1 << (7 - (bit_idx % 8)))) >> (7 - (bit_idx % 8))


def bits_from_bytearray(arr: bytearray, start_idx: int, num_bits: int, order: str) -> int:
    out = 0
    for bit_idx in range(start_idx, start_idx + num_bits):
        out = (out << 1) | bit_from_bytearray(arr, bit_idx, order)
    return out


def get_signed_value(bits: int, num_bits: int) -> int:
    # EXTEND
    # Table F.1
    # Figure F.12

    # if output should be positive, output our bits
    # if output should be negative, output the maximum negative value plus our bits
    if bits < 2**(num_bits - 1):  # Check if bits is less than the middle value
        min_val = (-1 << num_bits) + 1
        return min_val + bits

    return bits


def get_next_huffman_value(data: bytearray, data_pos: int, huff_table: Dict[Tuple[int, int], int],
                           debug_print=False) -> Tuple[int, int]:
    # DECODE
    encoded_bits = bit_from_bytearray(data, data_pos, "big")
    start_bit = data_pos
    curr_pos = data_pos + 1

    while (encoded_bits, curr_pos - start_bit) not in huff_table:
        encoded_bits = (encoded_bits << 1) | bit_from_bytearray(scan_data, curr_pos, "big")
        curr_pos += 1

    num_bits = curr_pos - start_bit

    if debug_print:
        print(f"encoded: {encoded_bits:0{num_bits}b}, length: {num_bits}")
    return huff_table[(encoded_bits, num_bits)], num_bits


def ycbcr_to_rgb(lum, chrom_blue, chrom_red) -> Tuple[int, int, int]:
    # https://www.impulseadventure.com/photo/jpeg-color-space.html
    red = chrom_red * (2 - 2 * 0.299) + lum
    blue = chrom_blue * (2 - 2 * 0.114) + lum
    green = (lum - 0.114 * blue - 0.299 * red) / 0.587
    return min(255, max(0, round(red + 128))), min(255, max(0, round(green + 128))), min(255, max(0, round(blue + 128)))


# Define JPEG segment names according to offset from 0xFFC0
# Table B.1
seg_names = ["SOF0 - Baseline DCT; Huffman", "SOF1 - Extended sequential DCT; Huffman",
             "SOF2 - Progressive DCT; Huffman", "SOF3 - Lossless (sequential); Huffman",
             "DHT - Define Huffman table(s)", "SOF5 - Differential sequential DCT; Huffman",
             "SOF6 - Differential progressive DCT; Huffman", "SOF7 - Differential lossless (sequential); Huffman",
             "JPG - Reserved for JPEG extensions", "SOF9 - Extended sequential DCT; Arithmetic",
             "SOF10 - Progressive DCT; Arithmetic", "SOF11 - Lossless (sequential); Arithmetic",
             "DAC - Define arithmetic coding conditioning(s)", "SOF13 - Differential sequential DCT; Arithmetic",
             "SOF14 - Differential progressive DCT; Arithmetic",
             "SOF15 - Differential lossless (sequential); Arithmetic",
             "RST0 - Restart with modulo 8 count 0", "RST0 - Restart with modulo 8 count 1",
             "RST0 - Restart with modulo 8 count 2", "RST0 - Restart with modulo 8 count 3",
             "RST0 - Restart with modulo 8 count 4", "RST0 - Restart with modulo 8 count 5",
             "RST0 - Restart with modulo 8 count 6", "RST0 - Restart with modulo 8 count 7",
             "SOI - Start of image", "EOI - End of image", "SOS - Start of scan", "DQT - Define quantization table(s)",
             "DNL - Define number of lines", "DRI - Define restart interval", "DHP - Define hierarchical progression",
             "EXP - Expand reference components", "JFIF header"]
seg_names.extend(["[Reserved: Application segments]"] * 16)
seg_names.extend(["[Reserved: JPEG extension]"] * 14)
seg_names.extend(["Comment", "[Invalid]"])

# Create necessary variables to store header information
huff_tables = {}
quant_tables = {}
sof: jpeg_headers.StartOfFrame

# Precalculate IDCT constants (A.3.3)
idct_lookup = []
for y in range(8):
    idct_row = []
    for x in range(8):
        uv_matrix = []
        for u in range(8):
            uv_row = []
            for v in range(8):
                cu = (1 / math.sqrt(2)) if u == 0 else 1
                cv = (1 / math.sqrt(2)) if v == 0 else 1
                uv_row.append(cu * cv *
                              math.cos(((2 * x + 1) * u * math.pi) / 16) * math.cos(((2 * y + 1) * v * math.pi) / 16))
            uv_matrix.append(uv_row)
        idct_row.append(uv_matrix)
    idct_lookup.append(idct_row)

with open(IMAGE_TO_OPEN, "rb") as f:
    block_id_bytes = f.read(2)
    while block_id_bytes:
        block_id = int.from_bytes(block_id_bytes, byteorder="big")
        pos = f.tell() - 2

        if block_id < 0xFFC0:
            print("Segment ID expected, not found.")
            break  # JPEG_SEG_ERR

        print(f"*** Marker: {seg_names[block_id - 0xFFC0]}, (0x{block_id:04X}) ***")
        print(f"OFFSET: {pos} (0x{pos:X})\n")

        if block_id == 0xFFD9:    # -------------------------------------- End of Image
            print(f"Decoding process took {time() - start_time} seconds")

            # Show decoded RGB array
            t = tkinter.Tk()
            t.resizable(False, False)
            t.title("decoded image")
            display = PhotoDisplay(t, sof.num_lines, sof.samples_per_line, image_rgb)
            t.mainloop()

            break
        elif block_id == 0xFFDA:  # -------------------------------------- Start of Scan
            sos = jpeg_headers.read_start_of_scan(f)
            assert sof is not None

            # Get maximum sampling factors and number of MCUs
            max_h_sampling_factor = 0
            max_v_sampling_factor = 0
            for component in sof.components:
                max_h_sampling_factor = max(max_h_sampling_factor, component.h_sampling_factor)
                max_v_sampling_factor = max(max_v_sampling_factor, component.v_sampling_factor)

            mcu_size_x = 8 * max_h_sampling_factor
            mcu_size_y = 8 * max_v_sampling_factor

            num_mcu_x = math.ceil(sof.samples_per_line / mcu_size_x)
            num_mcu_y = math.ceil(sof.num_lines / mcu_size_y)

            print(f"MCU Size: {mcu_size_x} x {mcu_size_y}")
            print(f"{num_mcu_x} MCU cols, {num_mcu_y} MCU rows\n")

            # Read data that's left
            start_pos = f.tell()
            f.seek(0, 2)
            end_pos = f.tell()

            f.seek(start_pos)
            scan_data = bytearray(f.read(end_pos - start_pos))

            # Scan for next marker and remove all stuff bytes in scan data
            marker_pos = None
            marker_pos_diff = 0
            for i in range(len(scan_data) - 2, 0, -1):
                # Remove stuff byte
                if scan_data[i:i + 2] == b'\xFF\x00':
                    scan_data.pop(i + 1)
                    marker_pos_diff += 1

                # Marker found if 0xFF exists without a stuff byte after
                elif scan_data[i:i + 2] > b'\xFF\x00':
                    marker_pos = i
                    marker_pos_diff = 0

            assert marker_pos is not None

            # Set scan data from start to next marker
            scan_data = scan_data[:marker_pos - marker_pos_diff]
            f.seek(start_pos + marker_pos)

            print("Scan data: (after bitstuff removed)")
            print("  " + "".join(f"{b:02x} " +
                                 ("\n  " if (idx + 1) % 36 == 0 else "") for idx, b in enumerate(scan_data[:720])))
            if len(scan_data) > 720:
                print("WARNING: Dump truncated.")
            print()
            # block_id_bytes = f.read(2)
            # continue

            # Initialize array storing final RGB values
            image_rgb = [[(0, 0, 0) for _ in range(sof.samples_per_line)] for _ in range(sof.num_lines)]

            # F.2.1.2
            # Figure E.9
            # Figure E.10
            curr_bit = 0

            # Initialize DC predictions for each component to zero
            predictions = [0 for _ in range(len(sos.components))]

            # Loop through all MCUs in image
            for mcu_row in range(num_mcu_y):
                print(f"Processing MCU row {mcu_row}")
                for mcu_col in range(num_mcu_x):
                    debug = False
                    if debug:
                        print(f"\nMCU {mcu_row}, {mcu_col}")

                    # List of MCUs (could contain luminance, chrominance blue, chrominance red or just one)
                    mcu_arr = []

                    # Decode each MCU
                    for component_idx, component in enumerate(sos.components):
                        frame_component = None
                        for c in sof.components:
                            if c.identifier == component.selector:
                                frame_component = c
                                break
                        assert frame_component is not None

                        quant_table = quant_tables[frame_component.quant_table_dest]
                        dc_huff_table = huff_tables[component.dc_table, 0]
                        ac_huff_table = huff_tables[component.ac_table, 1]

                        # Initialize 2D array for MCU
                        mcu = [[0 for _ in range(8 * frame_component.h_sampling_factor)]
                               for _ in range(8 * frame_component.v_sampling_factor)]

                        # Go through all data units in order specified by A.2.3
                        for data_unit_row in range(frame_component.v_sampling_factor):
                            for data_unit_col in range(frame_component.h_sampling_factor):
                                if debug:
                                    if frame_component.quant_table_dest == 0:
                                        print("Lum")
                                    else:
                                        print("Chr")

                                # Decode DC coefficient
                                dc_code, length = get_next_huffman_value(scan_data, curr_bit, dc_huff_table, debug)
                                curr_bit += length

                                # RECEIVE
                                # A.3.5, F.2.1.3.1
                                additional_bits = bits_from_bytearray(scan_data, curr_bit, dc_code, "big")
                                curr_bit += dc_code

                                diff = get_signed_value(additional_bits, dc_code)
                                abs_dc_value = predictions[component_idx] + diff

                                predictions[component_idx] = abs_dc_value

                                if debug:
                                    start_byte = (curr_bit - dc_code - length) // 8
                                    print(f"val: {diff}, coeff: 00=DC")
                                    print(f"val_bits: {additional_bits:0{dc_code}b}")
                                    print(f"data: 0x "
                                          f"{' '.join(hex(b)[2:] for b in scan_data[start_byte:start_byte + 4])}\n")

                                # Start decoding DCT matrix
                                dct_coeffs = [0 for _ in range(64)]
                                dct_coeffs[0] = abs_dc_value

                                # Decode AC coefficients, F.2.2.2
                                # Figure F.13
                                k = 0
                                while k != 63:
                                    k += 1

                                    # rs (8 bits)-> rrrrssss
                                    rs, length = get_next_huffman_value(scan_data, curr_bit, ac_huff_table, debug)
                                    curr_bit += length

                                    rrrr = rs >> 4  # Skip
                                    ssss = rs & 0b1111  # Coded length

                                    if ssss == 0:
                                        if rrrr == 15:
                                            k += 15
                                            continue
                                        else:
                                            if debug:
                                                print("EOB")
                                            break

                                    k += rrrr

                                    # Decode_ZZ
                                    # ZZ(k) = RECEIVE(ssss)
                                    additional_bits = bits_from_bytearray(scan_data, curr_bit, ssss, "big")
                                    curr_bit += ssss

                                    # ZZ(k) = EXTEND(ZZ(k), ssss)
                                    v = get_signed_value(additional_bits, ssss)

                                    if debug:
                                        start_byte = (curr_bit - ssss - length) / 8
                                        print(f"val: {v}, coeff: {k-rrrr:02d}..{k:02d}, skip: {rrrr}")
                                        print(f"val_bits: {additional_bits:0{ssss}b}")
                                        print(f"data: 0x "
                                              f"{' '.join(hex(b)[2:] for b in scan_data[start_byte:start_byte+4])}\n")

                                    dct_coeffs[k] = v

                                # Zig-zag reorder AC and DC coefficient list into DCT matrix
                                # Multiply each value by its corresponding value in the quantization table to dequantize
                                dct_matrix = [[0 for _ in range(8)] for _ in range(8)]
                                for i, coeff in enumerate(dct_coeffs):
                                    row, col = jpeg_headers.zigzag[i]
                                    dct_matrix[row][col] = coeff * quant_table[row][col]

                                # Perform IDCT (A.3.3)
                                for y in range(8):
                                    for x in range(8):
                                        val = 0
                                        for u in range(8):
                                            for v in range(8):
                                                val += idct_lookup[y][x][u][v] * dct_matrix[v][u]
                                        val /= 4

                                        # Assign value in MCU
                                        mcu[(data_unit_row * 8) + y][(data_unit_col * 8) + x] = val

                        # Expand MCU to maximum MCU size by duplicating values vertically or horizontally
                        horiz_multiplier = max_h_sampling_factor // frame_component.h_sampling_factor
                        vert_multiplier = max_v_sampling_factor // frame_component.v_sampling_factor

                        if vert_multiplier > 1 or horiz_multiplier > 1:
                            mcu = [[val for val in row for _ in range(horiz_multiplier)]
                                   for row in mcu for _ in range(vert_multiplier)]

                        # Append MCU to MCU array
                        mcu_arr.append(mcu)

                    # TODO: Handle images with one component
                    # Convert all Y, Cb, and Cr component values to RGB and store them in array
                    for i in range(mcu_size_y):
                        # Break if MCU goes past y bounds of image
                        if (mcu_row * mcu_size_y) + i >= sof.num_lines:
                            break

                        for j in range(mcu_size_x):
                            # Break if MCU goes past x bounds of image
                            if (mcu_col * mcu_size_x) + j >= sof.samples_per_line:
                                break

                            image_rgb[(mcu_row * mcu_size_y) + i][(mcu_col * mcu_size_x) + j] = \
                                ycbcr_to_rgb(mcu_arr[0][i][j], mcu_arr[1][i][j], mcu_arr[2][i][j])

        elif block_id == 0xFFD8:  # -------------------------------------- Start of Image
            pass
        elif block_id == 0xFFC4:  # -------------------------------------- Define Huffman Table
            for table in jpeg_headers.define_huffman_table(f):
                huff_tables[table.dest_id, table.table_class] = table.huff_data
        elif block_id == 0xFFDB:  # -------------------------------------- Define Quantization Table
            for table in jpeg_headers.define_quantization_table(f):
                quant_tables[table.dest_id] = table.table
        elif block_id == 0xFFC0 or block_id == 0xFFC1:  # ---------------- Start of Frame
            sof = jpeg_headers.read_start_of_frame(f)
        else:  # All other segments have length specified at the start, skip for now
            size = int.from_bytes(f.read(2), byteorder="big")
            f.seek(size - 2, 1)

        block_id_bytes = f.read(2)
