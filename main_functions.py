import os
import random
import logging


########################
####### ENCODING #######
########################
def encoding_picture(seed: int, num_bitlen: int, segment_len: int, rs_symbols_percent: float, correction_mode: str,  compression_mode: bool,
                     logging_filepath: str, logging_filename: str, logging_num_elements=100, logging_num_segments=10, rs_codec_value=32, min_codec_value=4,
                     min_segment_length=25, picture_path="Logo_long__binary_.png", picture_type="binary",
                     random_mask=True, show_image=True, verbose=True, do_xor=False):
    """
    This function encodes 1-bit images with the variants RS, RS_Spacer and DMR. The inserted RS symbols are calculated as a percentage of the total length of
    the section to be encoded. In addition, a logging file will be generated so that the individual steps of the encoding can be checked later.

    Args:
        seed: Start seed for the random mask.
        num_bitlen: Number of bits for byte translation.
        segment_len: Length of the DNA sequence a spacer should occur after.
        rs_symbols_percent: Percent of RS symbols from total fragment length. Ex: RS(255/223/32) 255 = 100%, 32 = 12.55% --> 0.1255 Input
        correction_mode: Determines the method with which the image is to be encoded. Selection is possible from RS, RS_spacer and DMR.
        compression_mode: Determines if the image should be compressed with packbits.
        logging_filepath: Defines the filepath for the logging file.
        logging_filename: Defines the filename of the logging file.
        logging_num_elements: Defines the number of elements which should be written in te logfile.
        logging_num_segments: Defines the number of segments which should be written in the logfile.
        picture_path: Specifies the path of the image to be encoded.
        picture_type: Defines what image type the selected image belongs to. Currently only binary images are supported.
        random_mask: Defines whether a random masked should be used to encode the picture.
        show_image: Determines whether the image should be displayed before translated.
        verbose: Sets how detailed the logging file should be recorded. True: very detailed. False: Warnings only

    Returns:
        Returns a list with the encoded DNA segments and a list with the binary translations.
    """

    import math
    import packbits
    import cv2 as cv
    import numpy as np
    import functions as func

    # Logging File
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:\n%(message)s\n',
        level=logging.DEBUG if verbose else logging.WARNING,
        datefmt='%d/%m/%Y %H:%M:%S',
        handlers=[
            logging.FileHandler(logging_filepath + logging_filename, mode='w'),
            logging.StreamHandler()
        ]
    )
    ROOT_DIR = os.path.abspath(os.curdir)
    im = cv.imread(ROOT_DIR + picture_path)
    h, w, c = im.shape

    logging.info("ENCODING INPUT: \n*PARAMETER:\n Seed: %s, bit length: %s, RS-percent: %s, segment size: %s, compression: %s, \n picture_mode: %s, "
                 "random_mask: %s, xor: %s \n*FILE: \n picture filepath: %s \n picture height: %s, picture width: %s, picture channel: %s"
                 % (seed, num_bitlen, rs_symbols_percent, segment_len, compression_mode, picture_type, random_mask, do_xor,
                    picture_path, h, w, c))

    # Daten einlesen
    if picture_type == "binary":
        picture = func.read_1bit_image(ROOT_DIR + picture_path, show_image)
        picture_str = "".join(str(item) for item in picture)
        logging.debug('Raw binary image (first %s elements out of %s):\n%s'
                      % (logging_num_elements, len(picture), "\n".join([picture_str[i:i + w] for i in range(0, logging_num_elements, w)])))

        packed_data = np.packbits(picture)
        logging.debug('Packed binary image (first %s elements out of %s):\n%s' % (logging_num_elements, len(packed_data), packed_data[0:logging_num_elements]))


    # Kompression durchführen
    compressed_data = np.array(list(packbits.encode(packed_data)), dtype="uint8") if compression_mode else packed_data
    if list(compressed_data) != list(packed_data):  # Logging Info
        logging.debug('Packbits compressed data (first %s elements out of %s):\n%s'
                      % (logging_num_elements, len(compressed_data), compressed_data[0:logging_num_elements]))

    if do_xor:
        random.seed(seed)
        xor_bits = list(func.bytes_to_bits(8, compressed_data)[0])
        xor_string = [random.randint(0, 1) for x in range(len(xor_bits))]
        xor = [int(xor_bits[item]) ^ xor_string[item] for item in range(len(xor_bits))]
        xor_bytes = func.bits_to_bytes(8, "".join(str(item) for item in xor))[0]
        compressed_data = xor_bytes
        random.seed()
        logging.debug("Generated Xor string (first %s values out of %s):\n%s"
                      % (logging_num_elements, len(xor_string), np.array(xor_string[0:logging_num_elements])))
        logging.debug("Picture xor generated xor string (first %s values out of %s):\n%s"
                      % (logging_num_elements, len(compressed_data), np.array(compressed_data[0:logging_num_elements])))

    # Bitlen bei größeren Galois Feld anpassen
    if num_bitlen != 8:
        compressed_data = func.bits_to_bytes(num_bitlen, func.bytes_to_bits(8, list(compressed_data))[0])[0]
        logging.debug('Data converted to %s bit integers (first %s elements out of %s):\n%s' % (
                       num_bitlen, logging_num_elements, len(compressed_data), np.array(compressed_data[0:logging_num_elements])))

    # random mask
    if random_mask:
        compressed_data = func.random_mask(compressed_data, seed, num_bitlen)
        logging.debug('Masked data (first %s elements out of %s):\n%s'
                      % (logging_num_elements, len(compressed_data), np.array(compressed_data[0:logging_num_elements])))
    else:
        logging.debug('No random Masked applied')


    if correction_mode == "No":

        logging.info('START NO ENCODING: TRANSLATING SEQUENCE')
        data_bits = func.bytes_to_bits(num_bitlen, compressed_data)[0]
        logging.debug("Bit sequence (first %s bits out of %s):\n%s" % (logging_num_elements, len(data_bits), data_bits[0:logging_num_elements]))


        original_dna = func.bits_to_dna("", data_bits)[0]
        logging.debug("DNA sequence translated with one-bit mapping (first %s bases out of %s):\n%s"
                      % (logging_num_elements, len(original_dna), original_dna[0:logging_num_elements]))

        # Output
        binary_segment_list = [data_bits]
        rs_message_list = [original_dna]

    elif correction_mode == "No_with_spacer":

        logging.info('START NO-WITH-SPACER ENCODING: TRANSLATING SEQUENCE')
        # Übersetzung in Bits, DNA
        binary_segment_list = [func.bytes_to_bits(8, list(segment))[0] for segment in compressed_data]
        logging.debug("Bit sequence (first %s segments out of %s):\n%s"
                      % (logging_num_segments, len(binary_segment_list), np.array(binary_segment_list[0:logging_num_segments])))


        rs_message_list = [func.bits_to_dna("", segment)[0] for segment in binary_segment_list]
        logging.debug("Input DNA segments (first %s segments out of %s):\n%s"
                      % (logging_num_segments, len(rs_message_list), np.array(rs_message_list[0:logging_num_segments])))


    # normaler üblicher RS Code mit Galois field von exp= 12 → 4095 damit ein Segment vorhanden in dem RS Symbole hinzugefügt werden
    # bei Logo_long_binary DMR: 4200 Segmente mit je 6 RS Symbolen --> = 25200 Symbole insgesamt
    # bei wetterdaten unterschiedlich? --> Wie wollen wir das da machen? %? oder codecvalue für jeden einzelnen einfügen?
    elif correction_mode == "RS":

        logging.info('START RS ENCODING: ADDING RS ERROR CORRECTING CODES AND TRANSLATING SEQUENCE')

        # add RS codes percent --> Achtung:  Prozentangabe von gesamtfragment also daten + rs angegeben !!!!
        # Bsp.: 16144 = 87.45 % Daten , 18461 = 100 % Daten + RS (87.45%+12,55%)  --> = length data/0,8745
        len_total = math.ceil(len(compressed_data)/(1-rs_symbols_percent))    # Länge Gesamtfragment
        symbol_num = math.ceil(rs_symbols_percent * len_total)   # 12,55 %*Gesamtlänge = 2317 RS Symbole
        data_encoded, max_errors, max_erasures = func.rs_encode_with_symbol_number([compressed_data], symbol_num, num_bitlen)
        logging.debug("This codec can correct up to %s errors and %s erasures simultaneously." % (max_errors, max_erasures))

        logging.debug("RS encoded data with 2^%s Galois Field (first %s elements out of %s):\n%s"
                      % (num_bitlen, logging_num_elements, len(data_encoded), np.array(data_encoded[0:logging_num_elements])))

        data_bits = func.bytes_to_bits(num_bitlen, data_encoded)[0]
        total_len_bits = len(data_bits)
        logging.debug("Bit sequence (first %s bits out of %s):\n%s" % (logging_num_elements, len(data_bits), data_bits[0:logging_num_elements]))


        original_dna = func.bits_to_dna("", data_bits)[0]
        total_len_dna = len(original_dna)
        logging.debug("DNA sequence translated with one-bit mapping (first %s bases out of %s):\n%s"
                      % (logging_num_elements, len(original_dna), original_dna[0:500]))

        # Output
        binary_segment_list = [data_bits]
        rs_message_list = [original_dna]


    elif correction_mode == "RS_spacer":   # z.B. Nutzung UV-Experiment!

        logging.info("START RS-SPACER ENCODING: ADDING RS ERROR CORRECTING CODES TO DNA-SEGMENTS AND TRANSLATING SEQUENCE")

        # beim prozentualen hinzufügen wird beim letzten Segment auch prozentual hinzugefügt --> besser immer die gleich anzahl --> einfacher beim decoden
        symbole_num = math.ceil(rs_symbols_percent * segment_len)  # Anzahl rs zeichen
        segment_in = int(segment_len - symbole_num)   # Länge der daten pro segment
        logging.info("Segmentgröße: " + str(segment_len) + ", Daten pro Segment: " + str(segment_in) + ", RS-Symbole pro Segment: " + str(symbole_num))
        data_segments = func.segment_string(compressed_data, segment_in)

        data_encoded, max_errors, max_erasures = func.rs_encode_with_symbol_number(data_segments, symbole_num, 8)
        logging.debug("RS encoded data with 2^%s Galois Field results in total length of %s (first %s segments out of %s):\n%s"
                      % (num_bitlen, sum(len(item) for item in data_encoded), logging_num_segments, len(data_encoded),
                         "\n".join(map(str, data_encoded[0:logging_num_segments]))))

        # Übersetzung in Bits, DNA
        binary_segment_list = [func.bytes_to_bits(8, list(segment))[0] for segment in data_encoded]
        logging.debug("Bit sequence with total length of %s (first %s segments out of %s):\n%s"
                      % (sum(len(item) for item in binary_segment_list), logging_num_segments, len(binary_segment_list),
                         np.array(binary_segment_list[0:logging_num_segments])))


        rs_message_list = [func.bits_to_dna("", segment)[0] for segment in binary_segment_list]
        logging.debug("Input DNA segments with total length of %s (first %s segments out of %s):\n%s"
                      % (sum(len(item) for item in rs_message_list), logging_num_segments, len(rs_message_list),
                         np.array(rs_message_list[0:logging_num_segments])))


    elif correction_mode == "DMR":

        import dmr_rs_coder as codec

        logging.info('TRANSLATING SEQUENCE USING DMR-RS, ADDING SPACERS')
        dmr_rsm_coder_main = codec.DMR_RS_Coder(rs_codec_value=rs_codec_value, min_codec_value=min_codec_value, min_segment_length=min_segment_length)

        rsm_segment_list, rs_symbol_size, payload_size = dmr_rsm_coder_main.encode_rsm(compressed_data)
        rsm_message_size = rs_symbol_size + payload_size
        logging.info("segments with: " + str(rsm_message_size * 8) + " length.")
        logging.debug("RS encoded data (first %s segments out of %s):\n%s"
                      % (logging_num_segments, len(rsm_segment_list), np.array(rsm_segment_list[0:logging_num_segments])))

        binary_segment_list = [func.bytes_to_bits(8, list(arr))[0] for arr in rsm_segment_list]
        logging.debug("Bit sequence (first %s segments out of %s. Total of %s elements.):\n%s"
                      % (logging_num_segments, len(binary_segment_list), sum(len(item) for item in binary_segment_list),
                         np.array(binary_segment_list[0:logging_num_segments])))

        dmr_message_list = [func.translate_binary_to_dmr_singular(binary_segment, segment_indice) for segment_indice, binary_segment in
                            enumerate(binary_segment_list)]
        logging.debug("Input DNA segments (first %s segments out of %s. Total of %s elements.):\n%s"
                      % (logging_num_segments, len(dmr_message_list), sum(len(item) for item in dmr_message_list),
                         np.array(dmr_message_list[0:logging_num_segments])))

        rs_message_list = dmr_message_list


    # Logger-Konfiguration für das zweite File
    logger = logging.getLogger()  # Aktuellen Logger holen
    for handler in logger.handlers:
        logger.removeHandler(handler)  # Vorhandene Handler entfernen

    return rs_message_list, binary_segment_list


def encoding_picture_segmented_packbits(seed: int, num_bitlen: int, segment_len: int, rs_symbols_percent: float, correction_mode: str, logging_filepath: str,
                                        logging_filename: str, logging_num_elements=100, logging_num_segments=10,
                                        picture_path="Logo_long__binary_.png", picture_type="binary", show_image=True, verbose=True):
    """
    This function encodes 1-bit images of the RS and DMR segmented packbit variant. Here, the image is read in and the strand of pixels is divided into segments
    of a certain size. These segments are then further encoded. The advantage of this is that the size of each segment is known, so that shifts due to packbit
    compression can be prevented.The inserted RS symbols are calculated as a percentage of the total length of the section to be encoded. In addition, a logging
    file will be generated so that the individual steps of the encoding can be checked later.

    Args:
        seed: Start seed for the random mask.
        num_bitlen: Number of bits for byte translation.
        segment_len: Length of the DNA sequence a spacer should occur after.
        rs_symbols_percent: Percent of RS symbols from total fragment length. Ex: RS(255/223/32) 255 = 100%, 32 = 12.55% --> 0.1255 Input
        logging_filepath: Defines the filepath for the logging file.
        logging_filename: Defines the filename of the logging file.
        logging_num_elements: Defines the number of elements which should be written in te logfile.
        logging_num_segments: Defines the number of segments which should be written in the logfile.
        correction_mode: Determines the method with which the image is to be encoded. Selection is possible from RS_segmented_packbits or DMR_segmented_packbits.
        picture_path: Specifies the path of the image to be encoded.
        picture_type: Defines what image type the selected image belongs to. Currently only binary images are supported.
        mapping_table: Choose between the mapping_table two_bit and no_random when not using DMR.
        show_image: Determines whether the image should be displayed before translating.
        verbose: Sets how detailed the logging file should be recorded. True: very detailed. False: Warnings only

    Returns:
        Returns a list with the encoded DNA segments and a list with the binary translations.
    """
    import math
    import packbits
    import cv2 as cv
    import numpy as np
    import functions as func

    # Logging File
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:\n%(message)s\n',
        level=logging.DEBUG if verbose else logging.WARNING,
        datefmt='%d/%m/%Y %H:%M:%S',
        handlers=[
            logging.FileHandler(logging_filepath + logging_filename, mode='w'),
            logging.StreamHandler()
        ]
    )

    ROOT_DIR = os.path.abspath(os.curdir)
    im = cv.imread(ROOT_DIR + picture_path)
    h, w, c = im.shape

    logging.info("ENCODING INPUT: \n*PARAMETER:\n Seed: %s, bit length: %s, RS-percent: %s, segment size: %s \n picture_mode: %s, "
                 "\n*FILE: \n picture filepath: %s \n picture height: %s, picture width: %s, picture channel: %s"
                 % (seed, num_bitlen, rs_symbols_percent, segment_len, picture_type, picture_path, h, w, c))

    # Daten einlesen
    if picture_type == "binary":
        picture = func.read_1bit_image(ROOT_DIR + picture_path, show_image)
        picture_str = "".join(str(item) for item in picture)
        logging.debug('Raw binary image (first %s elements out of %s):\n%s'
                      % (logging_num_elements, len(picture), "\n".join([picture_str[i:i + w] for i in range(0, logging_num_elements, w)])))

    # Segmentieren + Kompression
    segmented_picture = func.segment_string(picture, segment_len)
    packed_data = [np.packbits(segment) for segment in segmented_picture]
    packed_data_length = sum(len(item) for item in packed_data)
    logging.debug('Packed binary image (first %s segments out of %s. Total of %s elements.):\n%s'
                  % (logging_num_segments, len(packed_data), packed_data_length, packed_data[0:logging_num_segments]))


    compressed_data = [np.array(list(packbits.encode(segment)), dtype="uint8") for segment in packed_data]
    len_compressed_data = sum(len(compressed) for compressed in compressed_data)
    logging.info("Bild in %s lange Segmente geteilt und komprimiert. Insgesamt %s Segmente mit insgesamt %s Zeichen.\nFirst %s Segments: %s"
                 % (segment_len, len(compressed_data), len_compressed_data, logging_num_segments, compressed_data[0:logging_num_segments]))


    # Zufallsmaske hinzufügen
    masked = [func.random_mask(segment, seed + i, num_bitlen) for i, segment in enumerate(compressed_data)]
    len_masked_data = sum(len(mask) for mask in masked)
    logging.debug('Masked data (first %s segments out of %s. Total of %s elements.):\n%s'
                  % (logging_num_segments, len(masked), len_masked_data, masked[0:logging_num_segments]))


    if correction_mode == "RS_segmented_packbits":

        # add RS codes
        encoded = []
        for segment in masked:
            total_segment_len = math.ceil(len(segment)/(1 - rs_symbols_percent))   # Länge Gesamtfragment
            number_rs_symbols = math.ceil(rs_symbols_percent * total_segment_len)
            data_encoded, max_errors, max_erasures = func.rs_encode_with_symbol_number([segment], number_rs_symbols, num_bitlen)
            encoded.append(data_encoded)

        encoded_total_len = sum(len(encode) for encode in encoded)
        logging.debug("RS encoded data with 2^%s Galois Field (first %s segments out of %s. Total of %s elements.):\n%s"
                      % (num_bitlen, logging_num_segments, len(encoded), encoded_total_len, encoded[0:logging_num_segments]))


        # Übersetzung in Bits und DNA
        binary_segment_list = [func.bytes_to_bits(8, list(segment))[0] for segment in encoded]
        binary_seg_total_len = sum(len(b) for b in binary_segment_list)
        logging.debug("Bit sequence (first %s segments out of %s. Total of %s elements.):\n%s"
                      % (logging_num_segments, len(binary_segment_list), binary_seg_total_len,
                         np.array(binary_segment_list[0:logging_num_segments])))


        message_list = [func.bits_to_dna("", segment)[0] for segment in binary_segment_list]
        logging.debug("Input DNA segments (first %s segments out of %s. Total of %s elements.):\n%s"
                      % (logging_num_segments, len(message_list), sum(len(m) for m in message_list), np.array(message_list[0:logging_num_segments])))


    elif correction_mode == "DMR_segmented_packbits":
        import dmr_rs_coder as codec
        from reedsolo import RSCodec


        logging.info('TRANSLATING SEQUENCE USING DMR-RS')
        encoded = []
        rs_codec = rs_symbols_percent * 255
        for segment in masked:
            # falls prozentual hinzugefügt werden soll, können hiermit die Variablen bestimmt werden und in die Funktion eingegeben werden
            total_segment_len = int(math.ceil(len(segment) / (1 - rs_symbols_percent)))  # Länge Gesamtfragment
            number_RS_symbols = total_segment_len - len(segment)
            dmr_rsm_coder_main = codec.DMR_RS_Coder(rs_codec_value=rs_codec, min_codec_value=number_RS_symbols, min_segment_length=total_segment_len)

            rsm_segment_list, rs_symbol_size, payload_size = dmr_rsm_coder_main.encode_rsm(func)
            rsm_message_size = rs_symbol_size + payload_size
            encoded.append(rsm_segment_list[0])
        logging.debug("RS encoded data (first %s segments out of %s):\n%s" % (logging_num_segments, len(encoded), encoded[0:logging_num_segments]))


        binary_segment_list = [func.bytes_to_bits(8, list(arr))[0] for arr in encoded]
        logging.debug("Bit sequence (first %s segments out of %s. Total of %s elements.):\n%s"
                      % (logging_num_segments, len(binary_segment_list), sum(len(item) for item in binary_segment_list),
                         binary_segment_list[0:logging_num_segments]))


        message_list = [func.translate_binary_to_dmr_singular(binary_segment, segment_indice) for segment_indice, binary_segment in
                        enumerate(binary_segment_list)]
        logging.debug("Input DNA segments (first %s segments out of %s. Total of %s elements.):\n%s"
                      % (logging_num_segments, len(message_list), sum(len(item) for item in message_list), np.array(message_list[0:logging_num_segments])))


    logging.info("Encoding finished")

    # Logger-Konfiguration für das zweite File
    logger = logging.getLogger()  # Aktuellen Logger holen
    for handler in logger.handlers:
        logger.removeHandler(handler)  # Vorhandene Handler entfernen

    return message_list, binary_segment_list


########################
####### DECODING #######
########################
def picture_decode(dna_seq: [], seed: int, bitlen: int, rs_symbols_percent: float, segment_size: int, correction_mode: str, compression: bool,
                   picture_width: int, picture_height: int, save_filepath: str, save_filename: str, logging_filepath: str, logging_filename: str,
                   logging_num_elements=100, logging_num_segments=10, rs_codec_value=32, min_codec_value=4, min_segment_length=25,
                   picture_mode="binary", random_mask=True, verbose=True, do_xor=False, correction=True):
    """
    This function can be used to decode DNA sequences that were previously encoded using the encoding method: No, No_with_spacer, RS, RS_spacer or DMR in the
    picture_encode function. The inserted RS symbols are calculated as a percentage of the total length of the section to be encoded.
    In addition, a logging file will be generated so that the individual steps of the encoding can be checked later.

    Args:
        dna_seq: Input a list of dna sequences. The spacer need to be removed before decoding. Input the list of dna_sequences after spacer removal.
        seed: Start seed of the random mask.
        bitlen: Number of bits for byte translation.
        rs_symbols_percent: Percent of RS symbols from total fragment length. Ex: RS(255/223/32) 255 = 100%, 32 = 12.55% --> 0.1255 Input
        segment_size: Length of the DNA sequence a spacer should occur after.
        correction_mode: Determines the method to decode th image. Selection is possible from No, No_with_spacer, RS, RS_spacer or DMR.
        compression: Determines if the image was compressed with packbits.
        picture_width: Defines the width of the picture.
        picture_height: Defines the height of the picture.
        save_filepath: Defines the filepath for the output.
        save_filename: Defines the name for the output.
        logging_filepath: Defines the filepath for the logging file.
        logging_filename: Defines the filename of the logging file.
        logging_num_elements: Defines the number of elements which should be written in te logfile.
        logging_num_segments: Defines the number of segments which should be written in the logfile.
        picture_mode: Defines what image type the selected image belongs to. Currently only binary images are supported.
        random_mask: Defines whether the random masked should be used to decode.
        verbose: Sets how detailed the logging file should be recorded. True: very detailed. False: Warnings only

    Returns:
        Returns the corrected bit sequence and the picture. Furthermore, they are saved with the logging file in the specified path.

    """
    import math
    import packbits
    import numpy as np
    import functions as func
    from reedsolo import RSCodec

    # Logging File
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:\n%(message)s\n',
        level=logging.DEBUG if verbose else logging.WARNING,
        datefmt='%d/%m/%Y %H:%M:%S',
        handlers=[
            logging.FileHandler(logging_filepath + logging_filename, mode='w'),
            logging.StreamHandler()
        ]
    )

    logging.info("DECODING INPUT: \n*PARAMETER:\n Seed: %s, bit length: %s, RS-percent: %s, segment size: %s, compression: %s, \n picture_mode: %s, "
                 "picture-width: %s, picture-height: %s, \n random_mask: %s, xor: %s \n*FILE: \n filepath: %s \n filename: %s "
                 % (seed, bitlen, rs_symbols_percent, segment_size, compression, picture_mode, picture_width, picture_height, random_mask,
                    do_xor, save_filepath, save_filename))

    if correction_mode == "No":     # ohne RS

        logging.info('DECODING PICTURE WITH MODE NO')
        logging.debug("Original DNA sequence (first %s bases out of %s):\n%s"
                      % (logging_num_elements, len(dna_seq[0]), "\n".join([dna_seq[0][i:i + 75] for i in range(0, logging_num_elements, 75)])))

        # DNA to bits
        data_recovered = func.dna_to_bits(dna_seq[0])[0]
        logging.debug("Bit sequence (first %s bits out of %s):\n%s"
                      % (logging_num_elements, len(data_recovered), "\n".join([data_recovered[i:i + 75] for i in range(0, logging_num_elements, 75)])))

        # convert bits to bytes
        data_bytes = func.bits_to_bytes(bitlen, data_recovered)[0]
        logging.debug("Byte sequence (first %s elements out of %s):\n%s"
                      % (logging_num_elements, len(data_bytes), np.array(data_bytes[0:logging_num_elements])))
        binary_masked = data_bytes


    elif correction_mode == "No_with_spacer":

        logging.info('DECODING PICTURE WITH MODE NO-WITH-SPACER')
        # keine Spacer suchen --> rs_message_list als Eingabedatei --> erstmal einfacher!!!
        # Länge berechnen und ausgeben mit Hinweis + Spacer
        logging.debug("Original DNA sequence: total of %s elements(first %s segments out of %s):\n%s"
                      % (sum(len(item) for item in dna_seq), logging_num_segments, len(dna_seq), dna_seq[0:logging_num_segments]))

        # DNA to bits
        bits_list = [func.dna_to_bits(segment)[0] for segment in dna_seq]
        logging.debug("Recovered bit sequence (first %s segments out of %s):\n%s" % (logging_num_segments, len(bits_list), bits_list[0:logging_num_segments]))

        # convert bits to bytes
        bytes_list = [func.bits_to_bytes(bitlen, segment)[0] for segment in bits_list]
        logging.debug("Recovered byte sequences (first %s segments out of %s):\n%s"
                      % (logging_num_segments, len(bytes_list), np.array(bytes_list[0:logging_num_segments])))

        decoded_list = bytes_list
        # Länge der Decoded List Segmente anpassen → funktioniert nur für Logo-long nicht für Wetterdaten!
        # Bei den Wetterdaten würde das letzte Segment verlängert oder verkürzt werden --> wollen wir nicht

        for item in decoded_list:
            if len(item) > bitlen:  # Zeige Überschuss an Basen
                logging.error("Length of decoded list inaccurate. %s bases too many." % len(item[bitlen:]))
                del item[bitlen:]
            while len(item) < bitlen:
                logging.error("Length of decoded list inaccurate. %s bases too few." % (bitlen - len(item)))
                item.append(random.randint(0, 255))

        binary_masked = [item for sublist in decoded_list for item in sublist]
        logging.debug("Masked list (first %s bytes out of %s):\n%s"
                      % (logging_num_elements, len(binary_masked), np.array(binary_masked[0:logging_num_elements])))


    # normaler üblicher RS Code mit Galois field von exp= 12 → 4095 damit ein Segment vorhanden in dem RS Symbole hinzugefügt werden
    elif correction_mode == "RS":

        logging.info('DECODING PICTURE WITH MODE RS')
        logging.debug("Original DNA sequence (first %s bases out of %s):\n%s"
                      % (logging_num_elements, len(dna_seq[0]), "\n".join([dna_seq[0][i:i + 75] for i in range(0, logging_num_elements, 75)])))

        # DNA to bits
        data_recovered = func.dna_to_bits(dna_seq[0])[0]
        logging.debug("Bit sequence (first %s bits out of %s):\n%s"
                      % (logging_num_elements, len(data_recovered[0]), "\n".join([data_recovered[0][i:i + 75] for i in range(0, logging_num_elements, 75)])))

        # convert bits to bytes
        data_bytes = func.bits_to_bytes(bitlen, data_recovered)[0]
        logging.debug("Byte sequence (first %s elements out of %s):\n%s"
                      % (logging_num_elements, len(data_bytes), np.array(data_bytes[0:logging_num_elements])))

        rs_symbols = math.ceil(rs_symbols_percent * len(data_bytes))  # Anzahl rs zeichen vom Gesamtsegment --> hier jetzt Gesamtsegment
        if correction:
            # remove reed solomon codes
            rs_coder = RSCodec(rs_symbols, c_exp=bitlen)
            try:
                binary_masked = rs_coder.decode(data_bytes)[0]
                logging.debug("Byte sequence without RS Codes (first %s elements out of %s):\n%s" % (logging_num_elements,
                              len(binary_masked), np.array(binary_masked[0:logging_num_elements])))
            except:
                binary_masked = func.rs_force_uncorrect(rs_symbols, bitlen, data_bytes)[0]
                logging.debug("Too many errors for RS to correct, removing the appropriate number of symbols forcefully.\nByte sequence without RS Codes "
                              "(first %s elements out of %s):\n%s"
                              % (logging_num_elements, len(binary_masked), binary_masked[0:logging_num_elements]))

        else:    # keine Korrektur gewünscht
            binary_masked = func.rs_force_uncorrect(rs_symbols, bitlen, data_bytes)[0]
            logging.debug("No Correction chosen, removing the appropriate number of symbols forcefully.\nByte sequence without RS Codes "
                          "(first %s elements out of %s):\n%s"
                          % (logging_num_elements, len(binary_masked), np.array(binary_masked[0:logging_num_elements])))


    # abgewandelter RS für bessere DMR Vergleichbarkeit --> pro Segment sollen die RS Zeichen hinzugefügt werden -->
    # Segmente werden durch Spacer getrennt
    elif correction_mode == "RS_spacer":

        logging.info('DECODING PICTURE WITH MODE RS-SPACER')

        # keine Spacer suchen --> rs_message_list als Eingabedatei --> erstmal einfacher!!!
        # Länge berechnen und ausgeben mit Hinweis + Spacer
        dna_length = sum(len(item) for item in dna_seq)
        logging.debug("Original DNA sequence: (%s, first %s segments out of %s):\n%s"
                      % (dna_length, logging_num_segments, len(dna_seq), dna_seq[0:logging_num_segments]))

        # DNA to bits
        bits_list = [func.dna_to_bits(segment)[0] for segment in dna_seq]

        logging.debug("Recovered bit sequence (first %s segments out of %s):\n%s" % (logging_num_segments, len(bits_list), bits_list[0:logging_num_segments]))

        # convert bits to bytes
        bytes_list = [func.bits_to_bytes(bitlen, segment)[0] for segment in bits_list]
        logging.debug("Recovered byte sequences (first %s segments out of %s):\n%s"
                      % (logging_num_segments, len(bytes_list), bytes_list[0:logging_num_segments]))

        # remove reed solomon codes
        decoded_list = []
        rs_symbols = math.ceil(rs_symbols_percent * (segment_size / 8))  # Anzahl rs zeichen  # segmentsize als DNA Segment angegeben --> Byte benötigt / 8
        if correction:
            rs_coder = RSCodec(rs_symbols)
            for bytes_obj in bytes_list:
                try:
                    decoded_list.append(list(rs_coder.decode(bytes_obj)[0]))
                except:
                    forced = func.rs_force_uncorrect(rs_symbols, 8, bytes_obj)[0]
                    decoded_list.append(forced)
                    if len(forced) >= 10:
                        logging.error("Too many errors for RS to correct, removing the appropriate number of symbols forcefully."
                                      "\nByte sequence without RS Codes (first %s elements out of %s):%s"
                                      % (logging_num_segments, len(forced), forced[0:logging_num_segments]))

        else:   # keine Korrektur gewünscht
            for bytes_obj in bytes_list:
                forced = func.rs_force_uncorrect(rs_symbols, 8, bytes_obj)[0]
                decoded_list.append(forced)

        logging.debug("RS code removed (first %s segments out of %s):\n%s"
                      % (logging_num_segments, len(decoded_list), decoded_list[0:logging_num_segments]))

        # Länge der Decoded List Segmente anpassen → funktioniert nur für Logo-long nicht für Wetterdaten!
        # Bei den Wetterdaten würde das letzte Segment verlängert oder verkürzt werden --> wollen wir nicht
        len_data = int(segment_size/8 - rs_symbols)
        for item in decoded_list[0:-1]:
            if len(item) > len_data:   # Zeige Überschuss an Basen
                logging.error("Length of decoded list inaccurate. %s bases too many." % len(item[len_data:]))
                del item[bitlen:]
            while len(item) < len_data:
                logging.error("Length of decoded list inaccurate. %s bases too few." % (len_data - len(item)))
                item.append(random.randint(0, 255))

        binary_masked = [item for sublist in decoded_list for item in sublist]
        logging.debug(
            "Masked list (first %s bytes out of %s):\n%s" % (logging_num_elements, len(binary_masked), np.array(binary_masked[0:logging_num_elements])))


    elif correction_mode == "DMR":

        logging.info('DECODING PICTURE WITH MODE DMR')
        import dmr_rs_coder as codec
        from reedsolo import RSCodec

        if correction:
            dmr_rsm_coder_main = codec.DMR_RS_Coder(rs_codec_value=rs_codec_value, min_codec_value=min_codec_value, min_segment_length=min_segment_length)
            new_codec_size, payload_size = dmr_rsm_coder_main.recalculate_codec()
            errornous_segments, decodable_segments = dmr_rsm_coder_main.initial_scan_correction(dna_seq)
            errornous_segments_2, decodable_segments_2 = dmr_rsm_coder_main.after_scan_correction_try_Jess(errornous_segments)
            logging.debug("%s segment(s) decoded in initial scan. %s segment(s) decoded in after scan" % (len(decodable_segments), len(decodable_segments_2)))

            if len(errornous_segments_2) != 0:
                logging.error("%s segment(s) cannot be decoded by DMR:\n%s" % (len(errornous_segments_2), errornous_segments_2))

            all_segments = [entry for entry in decodable_segments]

            for entry in decodable_segments_2:
                all_segments += [entry]

            # Ansonsten würden die Segmente fehlen und müssten zufällig ersetzt werden → dann würde man mehr Fehler
            # erzeugen oder könnte seine Daten gar nicht ausgeben
            for entry in errornous_segments_2:
                all_segments += [entry]

            # Listen müssen sortiert werden damit die richtige Reihenfolge wieder da ist!
            all_segments.sort(key=lambda x: x[1])

            bits_list = [func.enhanced_translate_dna_to_binary_singular(*segment) for segment in all_segments]  # equal to line 98
            logging.debug("Recovered bit sequences (first %s segments out of %s):\n%s" % (logging_num_segments, len(bits_list), bits_list[0:logging_num_segments]))

            # Teilweise waren die Einträge in der Bits List zu kurz oder zu lang → hiermit passt es
            for number, item in enumerate(bits_list):
                if len(item) > (new_codec_size + payload_size) * 8:
                    logging.error("Segment size inaccurate. %s bases too many." % abs(len(item) - (new_codec_size + payload_size) * 8))
                    bits_list[number] = item[0: (new_codec_size + payload_size) * 8]

                while len(item) < (new_codec_size + payload_size) * 8:
                    logging.error("Segment size inaccurate. %s bases too few." % abs(len(item) - (new_codec_size + payload_size) * 8))
                    item += str(random.randint(0, 1))
                    bits_list[number] = item

                    while len(item) < (new_codec_size + payload_size) * 8:
                        item += str(random.randint(0, 1))
                        bits_list[number] = item

            bytes_list = func.bits_to_bytes(8, *bits_list)
            logging.debug("Recovered byte sequences (first %s segments out of %s):\n%s"
                          % (logging_num_segments, len(bytes_list), np.array(bytes_list[0:logging_num_segments])))

            # remove reed solomon codes
            decoded_list = []
            for bytes_obj in bytes_list:
                try:
                    rs_coder = RSCodec(new_codec_size)
                    decoded_list.append(list(rs_coder.decode(bytes_obj)[0]))
                except:
                    forced = func.rs_force_uncorrect(new_codec_size, bitlen, bytes_obj)[0]
                    decoded_list.append(forced)
                    logging.error("Too many errors for RS to correct, removing the appropriate number of symbols forcefully.\nByte sequence without RS Codes "
                                  "(first %s elements out of %s):%s" % (logging_num_segments, len(forced), forced[0:logging_num_segments]))

            logging.debug("RS code removed (first %s segments out of %s):\n%s"
                          % (logging_num_segments, len(decoded_list), np.array(decoded_list[0:logging_num_segments])))

            # decoded_list items waren teilweise unterschiedlich lang --> einige hatten 9 Einträge statt 8
            for item in decoded_list:
                if len(item) > payload_size:
                    logging.error("Length of decoded list inaccurate. %s bases too many." % len(item[payload_size:]))
                    del item[payload_size:]
                while len(item) < payload_size:
                    logging.error("Length of decoded list inaccurate. %s bases too few." % (payload_size - len(item)))
                    item.append(random.randint(0, 255))

        else:   # keine Korrektur mit DMR
            bits_list = []
            for index, segment in enumerate(dna_seq):
                try:
                    translated = func.enhanced_translate_dna_to_binary_singular(segment, index)
                    bits_list.append(translated)
                except:
                    bits_list.append("0")
            logging.debug("Recovered bit sequences (first %s segments out of %s):\n%s" % (logging_num_segments, len(bits_list), bits_list[0:logging_num_segments]))

            bytes_list = func.bits_to_bytes(8, *bits_list)
            logging.debug("Recovered byte sequences (first %s segments out of %s):\n%s"
                          % (logging_num_segments, len(bytes_list), np.array(bytes_list[0:logging_num_segments])))


            # remove reed solomon codes
            decoded_list = []
            # segmentgröße * (1+errorpercent)
            total_num_bytes = math.ceil(segment_size * (1 + rs_symbols_percent))
            rs_symboles = math.ceil(total_num_bytes * rs_symbols_percent)
            for bytes_obj in bytes_list:
                forced = func.rs_force_uncorrect(rs_symboles, bitlen, bytes_obj)[0]
                decoded_list.append(forced)
            logging.debug("RS code removed (first %s segments out of %s):\n%s"
                          % (logging_num_segments, len(decoded_list), np.array(decoded_list[0:logging_num_segments])))


        binary_masked = [item for sublist in decoded_list for item in sublist]
        logging.debug("Masked list (first %s bytes out of %s):\n%s"
                      % (logging_num_elements, len(binary_masked), np.array(binary_masked[0:logging_num_elements])))

    if random_mask:
        unmasked = func.remove_mask(binary_masked, seed, bitlen)
        logging.debug("Unmasked list (first %s bytes out of %s):\n%s" % (logging_num_elements, len(unmasked), np.array(unmasked[0:logging_num_elements])))
    else:
        unmasked = binary_masked
        logging.debug("No random masked needed to be removed.")

    bitseq = func.bytes_to_bits(8, unmasked)[0]   # Korrigierte Bitsequenz

    # for data in to decode ende --> Abschluss
    logging.info("CONVERTING BYTES BACK TO PICTURE AND WRITING REFERENCE FILE")
    if bitlen != 8:
        unmasked = func.bits_to_bytes(8, func.bytes_to_bits(bitlen, list(unmasked))[0])[0]
        logging.debug('Data converted to 8 bit integers (first %s elements out of %s):\n%s'
                      % (logging_num_elements, len(unmasked), np.array(unmasked[0:logging_num_elements])))

    if do_xor:
        random.seed(seed)
        xor_bits = list(func.bytes_to_bits(8, unmasked)[0])
        xor_string = [random.randint(0, 1) for x in range(len(xor_bits))]
        xor = [int(xor_bits[item]) ^ xor_string[item] for item in range(len(xor_bits))]
        xor_bytes = func.bits_to_bytes(8, "".join(str(item) for item in xor))[0]
        unmasked = xor_bytes
        random.seed()
        logging.debug("Generated Xor string (first %s values out of %s):\n%s"
                      % (logging_num_elements, len(xor_string), np.array(xor_string[0:logging_num_elements])))
        logging.debug("Unpacked binary data xor generated xor string (first %s values out of %s):\n%s"
                      % (logging_num_elements, len(xor), np.array(xor[0:logging_num_elements])))

    if picture_mode == "binary":
        if compression:
            try:
                decompressed_data = np.array(list(packbits.decode(unmasked)), dtype="uint8")
            except:
                decompressed_data = np.array(list(packbits.decode(unmasked[0:-1])), dtype="uint8")
            logging.debug("Packbits decompressed data (first %s values out of %s):\n%s"
                          % (logging_num_elements, len(decompressed_data), np.array(decompressed_data[0:logging_num_elements])))
        else:
            decompressed_data = np.array(unmasked, dtype="uint8")

        unpacked_data = list(np.unpackbits(decompressed_data))
        unpacked_str = "".join(str(item) for item in unpacked_data)
        logging.debug("Unpacked binary data (first %s values out of %s):\n%s"
                      % (logging_num_elements, len(unpacked_data),
                         "\n".join([unpacked_str[i:i + picture_width] for i in range(0, logging_num_elements, picture_width)])))

        logging.debug("Width: %s\nHeight: %s\nTotal number of pixels: %s\nLength of input binary list: %s" % (
                      picture_width, picture_height, picture_width * picture_height, len(unpacked_data)))


        output = func.binary_to_1bit_image(unpacked_data, picture_width, picture_height, save_filepath, save_filename)


    # Logger-Konfiguration für das zweite File
    logger = logging.getLogger()  # Aktuellen Logger holen
    for handler in logger.handlers:
        logger.removeHandler(handler)  # Vorhandene Handler entfernen

    return bitseq, output   # Binärfile, Picture


def picture_decode_segmented_packbits(dna_seq: [], seed: int, bitlen: int, rs_symbols_percent: float, segment_size: int, correction_mode: str,
                                      compression: bool, picture_width: int, picture_height: int, save_filepath: str, save_filename: str, logging_filepath: str,
                                      logging_filename: str, logging_num_elements=100, logging_num_segments=10, picture_mode="binary", correction=True,
                                      verbose=True):
    """
    This function can be used to decode DNA sequences that were previously encoded using the encoding method: RS_segmented_packbits and DMR_segmented_packbits
    in the picture_encode function. The inserted RS symbols are calculated as a percentage of the total length of the section to be encoded.
    In addition, a logging file will be generated so that the individual steps of the encoding can be checked later.

    Args:
        dna_seq: Input a list of dna sequences. The spacer need to be removed before decoding. Input the list of dna_sequences after spacer removal.
        seed: Start seed of the random mask.
        bitlen: Number of bits for byte translation.
        rs_symbols_percent: Percent of RS symbols from total fragment length. Ex: RS(255/223/32) 255 = 100%, 32 = 12.55% --> 0.1255 Input
        segment_size: Length of the DNA sequence a spacer should occur after.
        correction_mode: Determines the method to decode th image. Selection is possible from No, No_with_spacer, RS, RS_spacer or DMR.
        compression: Determines if the image was compressed with packbits.
        picture_width: Defines the width of the picture.
        picture_height: Defines the height of the picture.
        save_filepath: Defines the filepath for the output.
        save_filename: Defines the name for the output.
        logging_filepath: Defines the filepath for the logging file.
        logging_filename: Defines the filename of the logging file.
        logging_num_elements: Defines the number of elements which should be written in te logfile.
        logging_num_segments: Defines the number of segments which should be written in the logfile.
        picture_mode: Defines what image type the selected image belongs to. Currently only binary images are supported.
        verbose: Sets how detailed the logging file should be recorded. True: very detailed. False: Warnings only

    Returns:
        Returns the corrected bit sequence and the picture. Furthermore, they are saved with the logging file in the specified path.
        corrected bit sequence, picture, DMR: errornous_segments_2, counter: [# RS corrected seg., # forced corrected seg., DMR: # initial scan corrected seg.,
        DMR: # after scan corrected seg., DMR: # not corrected seg.]
    """

    import math
    import packbits
    import numpy as np
    import functions as func
    from reedsolo import RSCodec

    # Logging File
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:\n%(message)s\n',
        level=logging.DEBUG if verbose else logging.WARNING,
        datefmt='%d/%m/%Y %H:%M:%S',
        handlers=[
            logging.FileHandler(logging_filepath + logging_filename, mode='w'),
            logging.StreamHandler()
        ]
    )
    # RS corrected, forced correction, DMR initial scan, DMR after scan, DMR not corrected
    counter = [0, 0, 0, 0, 0]
    not_corrected = []   # interessant für Bearbeitung der Level

    logging.info("DECODING INPUT: \n*PARAMETER:\n Seed: %s, bit length: %s, RS-percent: %s, segment size: %s, compression: %s, \n picture_mode: %s, "
                 "picture-width: %s, picture-height: %s, \n*FILE: \n filepath: %s \n filename: %s "
                 % (seed, bitlen, rs_symbols_percent, segment_size, compression, picture_mode, picture_width, picture_height,
                    save_filepath, save_filename))

    # Länge berechnen und ausgeben mit Hinweis + Spacer
    dna_length = sum(len(item) for item in dna_seq)
    logging.debug("Original DNA sequence: (%s, first %s segments out of %s):\n%s"
                  % (dna_length, logging_num_segments, len(dna_seq)-1, dna_seq[0:logging_num_segments]))


    if correction_mode == "RS_segmented_packbits":
        logging.info('DECODING PICTURE WITH MODE RS-SEGMENTED-PACKBITS')

        # DNA to bits
        bits_list = [func.dna_to_bits(segment)[0] for segment in dna_seq]
        logging.debug("Recovered bit sequence (first %s segments out of %s):\n%s" % (logging_num_segments, len(bits_list), bits_list[0:logging_num_segments]))


    # DMR
    elif correction_mode == "DMR_segmented_packbits":

        logging.info('DECODING PICTURE WITH MODE DMR-SEGMENTED-PACKBITS')
        import scripts.encoding.dmr_rs as codec
        import scripts.decoding.decode_dmr as dmr

        if correction:      # DMR correction
            decoded_list = []

            for i, dna_segment in enumerate(dna_seq):

                try:
                    # falls prozentual hinzugefügt werden soll, können hiermit die Variablen bestimmt werden und in die Funktion eingegeben werden
                    rs_codec_value = rs_symbols_percent * 255
                    total_num_bytes = int(len(dna_segment) / 8)
                    rs_symbols = math.ceil(total_num_bytes * rs_symbols_percent)
                    dmr_rsm_coder_main = codec.DMR_RS_Coder(rs_codec_value=rs_codec_value, min_codec_value=rs_symbols,
                                                            min_segment_length=total_num_bytes)
                    new_codec_size, payload_size = dmr_rsm_coder_main.recalculate_codec()

                    errornous_segments, decodable_segments = dmr_rsm_coder_main.initial_scan_correction([dna_segment], index=i)
                    if len(decodable_segments) == 1:
                        decoded_list.append(decodable_segments)
                        counter[2] += 1
                    else:
                        errornous_segments_2, decodable_segments_2 = dmr_rsm_coder_main.after_scan_correction_try_Jess(errornous_segments)

                        if len(decodable_segments_2) == 1:
                            decoded_list.append(decodable_segments_2)
                            counter[3] += 1
                        else:
                            decoded_list.append(errornous_segments_2)
                            not_corrected.append(errornous_segments_2[0][0])
                            counter[4] += 1
                except:
                    decoded_list.append([(dna_segment, i)])
                    not_corrected.append(errornous_segments_2[0][0])
                    counter[4] += 1

            # Entfernen der Listen aus decoded List
            all_segments = [entry[0] for entry in decoded_list]
            decoded_list = all_segments

            # Listen müssen sortiert werden damit die richtige Reihenfolge wieder da ist! --> falls except aufgetreten ist!
            decoded_list.sort(key=lambda x: x[1])

            bits_list = [dmr.enhanced_translate_dna_to_binary_singular(*segment) for segment in decoded_list]  # equal to line 98
            logging.debug("Recovered bit sequences (first %s segments out of %s):\n%s" % (logging_num_segments, len(bits_list), bits_list[0:logging_num_segments]))

        else:
            decoded_list = dna_seq
            bits_list = []
            for index, segment in enumerate(decoded_list):
                try:
                    translated = dmr.translate_dna_to_binary_singular(segment, index)   # die enhanced funktion korrigiert auch schon!
                    bits_list.append(translated)
                except:
                    bits_list.append("00")
            logging.debug("Recovered bit sequences (first %s segments out of %s):\n%s" % (logging_num_segments, len(bits_list), bits_list[0:logging_num_segments]))

    # für alle weitermachen
    bytes_list = func.bits_to_bytes(bitlen, *bits_list)
    logging.debug("Recovered byte sequences (first %s segments out of %s):\n%s" % (logging_num_segments, len(bytes_list), bytes_list[0:logging_num_segments]))

    if correction:
        # remove reed solomon codes
        decoded_list = []
        for bytes_obj in bytes_list:
            rs_symboles = math.ceil(len(bytes_obj) * rs_symbols_percent)
            try:
                rs_coder = RSCodec(rs_symboles)
                decoded_list.append(list(rs_coder.decode(bytes_obj)[0]))
                counter[0] += 1
            except:
                forced = func.rs_force_uncorrect(rs_symboles, 8, bytes_obj)[0]
                decoded_list.append(forced)
                counter[1] += 1
                if len(forced) >= 10:
                    logging.error("Too many errors for RS to correct, removing the appropriate number of symbols forcefully."
                                  "\nByte sequence without RS Codes (first %s elements out of %s):%s"
                                  % (logging_num_segments, len(forced), forced[0:logging_num_segments]))
    else:
        # remove reed solomon codes
        decoded_list = []
        for bytes_obj in bytes_list:
            rs_symboles = math.ceil(len(bytes_obj) * rs_symbols_percent)
            forced = func.rs_force_uncorrect(rs_symboles, 8, bytes_obj)[0]
            decoded_list.append(forced)

    logging.debug("RS code removed (first %s segments out of %s):\n%s" % (logging_num_segments, len(decoded_list), decoded_list[0:logging_num_segments]))


    # Zufallsmaske entfernen
    binary_list = [func.remove_mask(decoded_list[i], seed + i, bitlen) for i in range(len(decoded_list))]
    logging.debug("Recovered unmasked list (first %s segments out of %s):\n%s" % (logging_num_segments, len(binary_list), binary_list[0:logging_num_segments]))

    # Dekompression
    decompressed_list = []
    for segment in binary_list:
        try:
            decompressed_list.append(np.array(list(packbits.decode(segment)), dtype="uint8"))
        except:
            segment.append(0)
            decompressed_list.append(np.array(list(packbits.decode(segment)), dtype="uint8"))

    logging.debug("Decompressed segments (first %s segments out of %s):\n%s"
                  % (logging_num_segments, len(decompressed_list), decompressed_list[0:logging_num_segments]))

    # uint8 wegmachen
    unpacked_data = []
    for segment in decompressed_list:
        try:
            unpacked_data.append(list(np.unpackbits(segment)))
        except:
            together = np.array(segment, dtype=np.uint8)
            unpacked_data.append(np.unpackbits(together))


    # recovered auf richtige Größe überprüfen:
    corrected = []
    for item in unpacked_data:
        while len(item) < segment_size:
            item = np.append(item, np.random.choice([0, 1]))
        if len(item) > segment_size:
            item = item[0: segment_size]
            corrected.append(item)
        else:
            corrected.append(item)


    # Listen zusammenfassen
    together = [item for sublist in corrected for item in sublist]

    if picture_mode == "binary":

        logging.debug("Width: %s\nHeight: %s\nTotal number of pixels: %s\nLength of input binary list: %s" % (
            picture_width, picture_height, picture_width * picture_height, len(together)))

        logging.info("WRITING RECOVERED BINARY IMAGE")
        output = func.binary_to_1bit_image(together, picture_width, picture_height, save_filepath, save_filename)


    # Logger-Konfiguration für das zweite File
    logger = logging.getLogger()  # Aktuellen Logger holen
    for handler in logger.handlers:
        logger.removeHandler(handler)  # Vorhandene Handler entfernen

    return together, output, not_corrected, counter  # Binärfile, Picture
