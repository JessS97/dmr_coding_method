import os
import json
import math
import random
import cv2 as cv
import Levenshtein
import numpy as np
from PIL import Image
import scipy.stats as stats
from reedsolo import RSCodec


##############################
### Read and write images ####
##############################


def read_1bit_image(fp, show_image=True):
    """
    This function converts an image to a black and white image containing only one bit for color coding.
    In addition, the image is compressed with packbits.
    :param fp: Sets the file path where the image was placed. The selected image must be an RGB image.
    :param show_image: Decide whether the image should be displayed.
    :return: This function outputs a list of the values of the generated black and white image from the source image,
    which were compressed using packbits.
    """

    """1. Vorbereiten des Bildes: Umwandeln in Schwarz-Weiß und Hinzufügen in Liste."""
    thresh = thresh_image(fp, show_image)             # Wird in schwarz-weiß Bild umgewandelt
    image_np_array = np.array(thresh)     # converts PIL images into NumPy arrays
    image_list = image_np_array.tolist()  # NumPy Liste zu normaler Liste
    image_1_bit_list = []

    """2. Umwandeln des Schwarz Weiß Bildes in 1 Bit Image."""
    for row in image_list:
        for value in row:
            if value != 0:   # Wenn Value nicht gleich 0 dann = 1
                image_1_bit_list += [value // 255]  # // = floor division: rounds the result down to the nearest whole
            else:                                      # number [255//255] = 1, Restliche Werte = 0
                image_1_bit_list += [value]   # Value = 0 --> übernehmen

    return image_1_bit_list


def thresh_image(fp, show_image=True):
    """
    This function converts a color image to a grayscale image and then to a black and white image.
    The conversion to a black and white image is done with a threshold. Pixel values below a certain threshold are
    set to 0 and the rest to 255. The optimal threshold for the image is determined by the Otsu-algorithm.
    :param fp: Sets the file path where the image was placed. The selected image must be an RGB image.
    :param show_image: Decide whether the image should be displayed.
    :return: The function returns the black and white version of the input image.
    """

    """ 1. Read the image. """
    image = cv.imread(fp)

    """ 2. Convert an image from BGR to grayscale mode. """
    gray_image = cv.cvtColor(image, cv.COLOR_BGR2GRAY)

    """ 3. Convert a grayscale image to black and white using binary thresholding. """
    (thresh, BnW_image) = cv.threshold(gray_image, 20, 255, cv.THRESH_BINARY | cv.THRESH_OTSU)
    # (Graustufenbild, Threshold Value, max. Pixelwert, Threshold Funktion, Otsu algorithm: choose opt. threshold value)
    if show_image:
        cv.imshow("image", BnW_image)   # cv2.imshow() method is used to display an image in a window.
        cv.waitKey(0)   # Das Window des Bildes bleibt auf bis eine Taste gedrückt wird.

    return BnW_image


def binary_to_1bit_image(binary_list: list, width: int, height: int, filepath: str, filename: str):
    """
    This function plays back a 1-bit image from an entered binary list. It checks if the list contains entries and if
    the list size fits to the defined image size. List entries that are too much are deleted and not considered for the
    translation into the image. In addition, the image will be saved under the other specified path with the specified
    name.
    Args:
        binary_list: Defines the binary string to be translated into a 1-bit image.
        width: Defines the image width.
        height: Defines the image height.
        filepath: Defines the directory the image should be saved.
        filename: Defines the name under which the image should be saved.

    Returns: Returns the created Numpy_Array.
    """

    image_1_bit_list = []
    if len(binary_list) == 0:
        raise ValueError("The binary list entered does not contain any entries. Please change the entry.")

    # Prüfen, ob die Bildgröße passt. Falls zu groß soll der Überhang entfernt werden.
    pixel = width * height

    if pixel == 0:
        raise ValueError("The entered image size corresponds to 0x0 pixels. Please enter the correct dimensions.")

    if pixel < len(binary_list):
        data = binary_list[0: pixel]
    elif pixel > len(binary_list):
        data = binary_list
        for add in range(pixel-len(binary_list)):
            data.append(random.choice([0, 1]))
    else:
        data = binary_list

    # 0 und 1 in Graustufen (0 und 255) wiedergeben
    image_1_bit_list = [255 if element != 0 else 0 for element in data]

    # Liste mit Arrays der Bildgrößen erstellen
    picture_list = lambda image_1bit_list, wide: [image_1_bit_list[i:i + wide] for i in range(0, len(image_1_bit_list), wide)]
    final_picture_list = picture_list(image_1_bit_list, width)

    # in Numpy Array umwandeln
    image_np_array = np.array(final_picture_list)

    # Bild einlesen und speichern
    picture = Image.fromarray(image_np_array.astype('uint8'), mode='L')
    picture.save(f"{filepath}\\{filename}.png")

    return image_np_array



###########################
######## Bit, Byte ########
###########################


def bytes_to_bits(bit_number: int, *args: list):
    """
    Converts byte objects to bits in string form (string of 1's and 0's)

    Args:
        bit_number: how many bits to be converted into
        *args: objects that need to be converted into bit sequences

    Returns:
        For each object entered in *bytes_objects, the object as a string of 1's and 0's corresponding to the sequence
        of bits.
    """
    bits_list = [''.join([bin(item).replace("0b", "").zfill(bit_number) for item in arg]) for arg in args]

    return tuple(bits_list)


def bits_to_bytes(bit_number: int, *args: str):
    """
    Converts string of bits to byte object

    Args:
        bit_number: how many bits to be converted into
        *args: strings of bit sequences that need to be converted to bytes objects

    Returns:
        For each string entered in *bits_strings, the string as a bytes object
    """
    bytes_list = [[int(arg[i:i + bit_number], 2) for i in range(0, len(arg), bit_number)] for arg in args]

    return tuple(bytes_list)


###########################
###### Random Mask ########
###########################


def random_mask(value_array, seed, bitlen=8):
    """
    This function creates a random mask over an input array of values. This mask is used to vary the pixel values of
    a uniform image. Without this variation the image could not be stored in DNA, because no successful synthesis
    could take place.

    Args:
        value_array: Defines the values to be randomized. These must be entered as an array and can for example come from an image.
        seed: Defines the seed to be used by the random generator.
        bitlen: A positive integer that expresses the number of bits in a bit string.

    Returns:
        The function returns a list of randomized values.
    """

    random_gridded_value_list = []
    random.seed(seed)

    for value in value_array:
        random_gridded_value_list += [(value + random.randint(0, 2**bitlen-1)) % 2**bitlen]

    random.seed()

    return random_gridded_value_list


def remove_mask(decoded_list: list, seed: int, bitlen=8):
    """
    This function removes the random mask, which was used to vary the pixel values of a uniform image. Without this variation the image could not be stored in
    DNA, because no successful synthesis could take place.

    Args:
        decoded_list: A list with the corrected byte values.
        seed: Defines the seed to be used by the random generator.
        bitlen: A positive integer that expresses the number of bits in a bit string.

    Returns:
        The function returns a list with the unmasked values.
    """

    unmasked_value = []
    random.seed(seed)
    for value in decoded_list:
        fixed = value - random.randint(0, 2**bitlen-1)
        if fixed >= 0:
            unmasked_value.append(fixed)
        else:
            unmasked_value.append(fixed+2**bitlen)
    random.seed()

    return unmasked_value


###########################
###### Segmentation #######
###########################


def segment_string(message, segment_size):
    """
    This Function segments the given Message to Segments of n-length.
    It is mainly meant to be used for segmenting the DMR_String but can also be used to segment the binary Message.
    Args:
        message: The input can be a DNA strand as well as a binary sequence.
        segment_size: Defines the size of the segments into which the input sequence should be divided.
    Returns:
        Returns a list with the segmented data.
    """

    segment_quantity = math.ceil(len(message) / segment_size)       # Bestimmung der gesamten Segmentanzahl
    segment_list = []
    for i in range(segment_quantity):
        segment = message[i * segment_size: (i + 1) * segment_size]
        if len(segment) != 0:             # Falls keine Sequenz eingegeben wurde nicht zur Liste hinzufügen
            segment_list += [segment]     # Packt auch letzte kürzere Segmente hinzu

    return segment_list


#########################
###### Encode DMR #######
#########################


def translate_binary_to_dmr_singular(segment, segment_index):
    """
    This function translates a binary String to DNA on basis of the DMR-Scheme. The DMR-Scheme applies a certain
    pattern (Dynamic-Mapping-Rule) to the DNA that can be partially retraced for correction. In contrast to the
    translate_binary_to_DMR function, this translates a single fragment of the already segmented binary string.
    :param segment: Input of a segment of a segmented binary string.
    :param segment_index: Wie soll die Eingabe hier erfolgen?
    :return: Outputs the translated DNA strand.
    """

    """1. Aufrufen der Map Libraries und Vorbereitung """
    # Lade die JSON-Datei in ein Dictionary
    root_dir = os.path.abspath(os.curdir)
    with open(f"{root_dir}\\mapping_table_dmr.json") as f:
        data = json.load(f)

    map_library = data["map_library"]
    map_library_start_2_mere = data["initial_2mer"]

    # Variablen
    dmr_rs_dna_message = ""
    last_coded_2_mere = ""
    binary_pair = segment[0:2]            # Bestimmung des ersten 2-mers
    segment_modular = segment_index % 4

    # Für Einlesen der json Dateien. Konnte ohne nicht auf die Elemente zugreifen.
    if segment_modular == 0:
        segment_mod = "0"
    elif segment_modular == 1:
        segment_mod = "1"
    elif segment_modular == 2:
        segment_mod = "2"
    else:
        segment_mod = "3"

    # Wie soll das Segment eingegeben werden? Schon mod 4 oder verarbeitet oder normale Segmentnummer?
    # segment_count = segment_count % 4 hinzufügen?

    """2. Übersetzung """

    # <editor-fold desc="2.1 Übersetzen des Start 2-meres mit der Mappingtabelle 1 und dem jeweiligen Segment mod 4.">
    if binary_pair == "00":     # zuvor bestimmtes Start 2-mere
        two_mere = map_library_start_2_mere[segment_mod][0]       # von Mapping Tabelle übernehmen
        dmr_rs_dna_message += two_mere                              # zur Sequenz hinzufügen
        last_coded_2_mere = two_mere                                # 2- mere merken

    elif binary_pair == "01":
        two_mere = map_library_start_2_mere[segment_mod][1]
        dmr_rs_dna_message += two_mere
        last_coded_2_mere = two_mere

    elif binary_pair == "10":
        two_mere = map_library_start_2_mere[segment_mod][2]
        dmr_rs_dna_message += two_mere
        last_coded_2_mere = two_mere

    elif binary_pair == "11":
        two_mere = map_library_start_2_mere[segment_mod][3]
        dmr_rs_dna_message += two_mere
        last_coded_2_mere = two_mere
    # </editor-fold>

    # <editor-fold desc="2.2. Übersetzen der restlichen Basen nach Mapping Tabelle 2 unter Einbezug der Vorherigen.">

    for j in range(2, len(segment), 2):
        binary_pair = segment[j: j + 2]

        if binary_pair == "00":
            two_mere = map_library[last_coded_2_mere][0]
            dmr_rs_dna_message += two_mere
            last_coded_2_mere = two_mere

        elif binary_pair == "01":
            two_mere = map_library[last_coded_2_mere][1]
            dmr_rs_dna_message += two_mere
            last_coded_2_mere = two_mere

        elif binary_pair == "10":
            two_mere = map_library[last_coded_2_mere][2]
            dmr_rs_dna_message += two_mere
            last_coded_2_mere = two_mere

        elif binary_pair == "11":
            two_mere = map_library[last_coded_2_mere][3]
            dmr_rs_dna_message += two_mere
            last_coded_2_mere = two_mere

    # </editor-fold>

    return dmr_rs_dna_message


#########################
###### Decode DMR #######
#########################


def translate_dna_to_binary_singular(segment, index):
    """
    This function translates DNA, which was previously translated from a binary sequence by the DMR scheme,
    back into a binary sequence. The back translation is also done with the help of the DMR scheme. In addition,
    the actual 2-mers is compared with the according to the rule possible 2-mers.
    :param segment: Input of the segment of a DNA sequence.
    :param index: Definition of the segment number of the entered DNA strand.
    :return: Returns a binary sequence.
    """

    """1. Aufrufen der Map Libraries und Vorbereitung """
    segment_modular = index % 4  # Bestimmung der Segmentnummer für das DMR Schema
    # Für Einlesen der json Dateien. Konnte ohne nicht auf die Elemente zugreifen.
    if segment_modular == 0:
        segment_count = "0"
    elif segment_modular == 1:
        segment_count = "1"
    elif segment_modular == 2:
        segment_count = "2"
    else:
        segment_count = "3"

    # Lade die JSON-Datei in ein Dictionary
    root_dir = os.path.abspath(os.curdir)
    with open(f"{root_dir}\\mapping_table_dmr.json") as f:
        data = json.load(f)
    # Zugriff auf den jeweiligen Teil des Dictionaries
    map_library = data["map_library"]
    map_library_start_2_mere = data["initial_2mer"]

    binary_segment = ""
    start_2_mere = segment[0:2]  # Beschreibt das erste 2-mere des Segments
    position_in_rule = None  # Beschreibt die Position in jeweiligen Zeile der Regel, Position 1: [00], 2. [01],
    last_2_mere = None  # Merkt sich das letzte 2-mere

    """2. Übersetzung """

    # <editor-fold desc="2.1 Übersetzung des 1. 2-meres nach der 1. Mapping Tabelle.">
    # Positionsbestimmung der jeweiligen Zeile in der Tabelle
    for i, rule_2_mere in enumerate(map_library_start_2_mere[segment_count]):  # enumerate gibt Count und Value zurück
        # Ausgabe Bsp.: Segment count = 0: 1. Count [00] = AA, 2. Count [01] = CC, 3. [10] = GG, 4. [11] = TT

        if rule_2_mere == start_2_mere:  # geht die Sachen einzeln durch und vergleicht diese. Wenn identisch, dann..
            position_in_rule = i  # Dies ist die Position in der Regel
            last_2_mere = rule_2_mere  # 2-mere als letztes 2-mere merken
            break  # nicht weiter durchgehen

    # Übersetzen des 2-meres anhand der bestimmten Position
    if position_in_rule == 0:
        binary_segment += "00"

    elif position_in_rule == 1:
        binary_segment += "01"

    elif position_in_rule == 2:
        binary_segment += "10"

    elif position_in_rule == 3:
        binary_segment += "11"

    current_2_mere = None  # Kann das raus?

    # </editor-fold>

    # <editor-fold desc="2.1 Übersetzung der restlichen Basen mithilfe der 2. Mapping Tabelle.">

    for i in range(2, len(segment), 2):
        current_2_mere = segment[i:i + 2]

        for j, rule_2_mere in enumerate(map_library[last_2_mere]):
            if rule_2_mere == current_2_mere:
                position_in_rule = j
                last_2_mere = rule_2_mere
                break

        if position_in_rule == 0:
            binary_segment += "00"

        elif position_in_rule == 1:
            binary_segment += "01"

        elif position_in_rule == 2:
            binary_segment += "10"

        elif position_in_rule == 3:
            binary_segment += "11"
    # print(binary_segment)
    # </editor-fold>

    return binary_segment


def enhanced_translate_dna_to_binary_singular(segment, index, verbose=True):
    """
    This function translates DNA, which was previously translated from a binary sequence by the DMR scheme,
    back into a binary sequence. This enhanced version now tries out the normal translation with the mapping scheme
    first. If no translation would be possible with it, it is tried to go through the variations of the last_2_meres
    whether the next 2-mere is found as rule in one of the variations. If a variation contains the next 2-mere as a
    rule, it is assumed that the variation is the actual current-2-mere and the translation is performed with it. If
    no next2-mere is found as a variation of the last 2-meres (e.g., in case of 2 wrong 2 mers in a row) then 6 random
    0s or 1s are added and the last-2-mere is set to the 2-mere 3 iterations removed. Likewise, 2 iterations are skipped,
     so that after the random insertion at the correct place can be translated further.
    :param segment: Input of the segment of a DNA sequence.
    :param index: Definition of the segment number of the entered DNA strand.
    :return: Returns a binary sequence.
    """

    """1. Aufrufen der Map Libraries und Vorbereitung """
    segment_modular = index % 4  # Bestimmung der Segmentnummer für das DMR Schema
    # Für Einlesen der json Dateien. Konnte ohne nicht auf die Elemente zugreifen.
    if segment_modular == 0:
        segment_count = "0"
    elif segment_modular == 1:
        segment_count = "1"
    elif segment_modular == 2:
        segment_count = "2"
    else:
        segment_count = "3"

    # Lade die JSON-Datei in ein Dictionary
    root_dir = os.path.abspath(os.curdir)
    with open(f"{root_dir}\\mapping_table_dmr.json") as f:
        data = json.load(f)
    # Zugriff auf den jeweiligen Teil des Dictionaries
    map_library = data["map_library"]
    map_library_start_2_mere = data["initial_2mer"]

    binary_segment = ""
    start_2_mere = segment[0:2]  # Beschreibt das erste 2-mere des Segments
    position_in_rule = None  # Beschreibt die Position in jeweiligen Zeile der Regel, Position 1: [00], 2. [01],
    last_2_mere = None  # Merkt sich das letzte 2-mere

    """2. Übersetzung """

    # <editor-fold desc="2.1 Übersetzung des 1. 2-meres nach der 1. Mapping Tabelle.">
    # Positionsbestimmung der jeweiligen Zeile in der Tabelle
    for i, rule_2_mere in enumerate(map_library_start_2_mere[segment_count]):  # enumerate gibt Count und Value zurück
        # Ausgabe Bsp.: Segment count = 0: 1. Count [00] = AA, 2. Count [01] = CC, 3. [10] = GG, 4. [11] = TT

        if rule_2_mere == start_2_mere:  # geht die Sachen einzeln durch und vergleicht diese. Wenn identisch, dann..
            position_in_rule = i  # Dies ist die Position in der Regel
            last_2_mere = rule_2_mere  # 2-mere als letztes 2-mere merken
            break  # nicht weiter durchgehen

    # Übersetzen des 2-meres anhand der bestimmten Position
    if position_in_rule == 0:
        binary_segment += "00"

    elif position_in_rule == 1:
        binary_segment += "01"

    elif position_in_rule == 2:
        binary_segment += "10"

    elif position_in_rule == 3:
        binary_segment += "11"

    current_2_mere = None  # Kann das raus?

    # </editor-fold>

    # <editor-fold desc="2.2 Übersetzung der restlichen Basen mithilfe der 2. Mapping Tabelle.">
    not_in_scheme = []
    continue_again = False
    continue_again_again = False
    for i in range(2, len(segment), 2):

        if continue_again:
            continue_again = False
            continue
        if continue_again_again:
            continue_again_again = False
            continue

        current_2_mere = segment[i:i + 2]
        wanted_length = i + 2

        # 2.2.1 normale Übersetzung ausprobieren.
        try:

            for j, rule_2_mere in enumerate(map_library[last_2_mere]):
                if rule_2_mere == current_2_mere:
                    position_in_rule = j
                    last_2_mere = rule_2_mere
                    break
                else:
                    position_in_rule = 4

            if position_in_rule == 0:
                binary_segment += "00"

            elif position_in_rule == 1:
                binary_segment += "01"

            elif position_in_rule == 2:
                binary_segment += "10"

            elif position_in_rule == 3:
                binary_segment += "11"

            elif position_in_rule == 4:
                raise ValueError("End try")


        except:  # normale Übersetzung funktioniert nicht

            # 2.2.2 Variation des Mapping Schemas probieren
            try:

                # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in map_library[last_2_mere]:
                    for j, rule_2_mere in enumerate(map_library[variation]):
                        next_2_mere = segment[i + 2:i + 4]

                        # Rule in Variation stimmt überein → rule_2_mere als current_2_mere annehmen
                        if rule_2_mere == next_2_mere:  # variation ist der eigentliche 2-mere
                            current_2_mere = variation
                            for j, rule_2_mere in enumerate(map_library[last_2_mere]):
                                if rule_2_mere == current_2_mere:
                                    position_in_rule = j
                                    last_2_mere = rule_2_mere
                                    break

                            if position_in_rule == 0:
                                binary_segment += "00"
                                break

                            elif position_in_rule == 1:
                                binary_segment += "01"
                                break

                            elif position_in_rule == 2:
                                binary_segment += "10"
                                break

                            elif position_in_rule == 3:
                                binary_segment += "11"
                                break

                            break

                    if len(binary_segment) == wanted_length:
                        break

                # Wenn nicht gepasst hat, wurde nichts hinzugefügt --> weiter mit except
                if len(binary_segment) < wanted_length:
                    raise ValueError("End try2")

            # Keine Variation durch Mapping Schema funktioniert:
            except:
                # zwei Fehler hintereinander --> nächstes 2-mere konnte aus der Mapping Tabelle nicht gefunden werden
                # Bsp: richtige Sequenz: TT TT TT AC AC
                #  falsche Sequenz:  TT TA TT AC AC
                #                        00 XX XX ab hier weiter
                if len(binary_segment) >= 8:
                    binary_segment += binary_segment[-8: -2]
                    not_in_scheme.append(i + 1)
                    not_in_scheme.append(i + 3)
                    not_in_scheme.append(i + 5)
                elif len(binary_segment) < 8:
                    for loop in range(6):
                        binary_segment += str(random.randint(0, 1))
                    not_in_scheme.append(i + 1)
                    not_in_scheme.append(i + 3)
                    not_in_scheme.append(i + 5)

                last_2_mere = segment[i + 4: i + 6]
                continue_again = True  # Überspringt die nächste for Iteration, da bereits hinzugefügt.
                continue_again_again = True

    # </editor-fold>

    return binary_segment



def trans_map_lib_one_2_mere(last_2_mere: str, current_2_mere: str):
    # Lade die JSON-Datei in ein Dictionary
    root_dir = os.path.abspath(os.curdir)
    with open(f"{root_dir}\\mapping_table_dmr.json") as f:
        data = json.load(f)
    # Zugriff auf den jeweiligen Teil des Dictionaries
    map_library = data["map_library"]

    for j, rule_2_mere in enumerate(map_library[last_2_mere]):
        if rule_2_mere == current_2_mere:
            position_in_rule = j
            break

    # Übersetzen des 2-meres anhand der bestimmten Position
    if position_in_rule == 0:
        binary_segment = "00"
    elif position_in_rule == 1:
        binary_segment = "01"
    elif position_in_rule == 2:
        binary_segment = "10"
    elif position_in_rule == 3:
        binary_segment = "11"

    return binary_segment


def dmr_check_mapping(sequence, segment_number=0, start=True):
    """
    This function checks whether a DNA sequence fits into the DMR scheme and could have been encoded with it.
    :param sequence: The DNA sequence to be compared should be 10 bases long and entered as a string.
    :param segment_number: Defines which segment we are in the DMR scheme for mapping.
    :return:  If the DMR scheme matches, true is returned. Otherwise false.
    """

    # Lade die JSON-Datei in ein Dictionary
    root_dir = os.path.abspath(os.curdir)
    with open(f"{root_dir}\\mapping_table_dmr.json") as f:
        data = json.load(f)
    # Zugriff auf den jeweiligen Teil des Dictionaries
    map_library = data["map_library"]
    start_library = data["initial_2mer"]

    decision_support = []  # Entscheidungshilfe zum Merken, ob die Einträge gleich waren

    if start:
        counter_2mere = 1  # Zähler welcher 2-mere gerade betrachtet wird
        """Start 2-mere muss in erster Tabelle vorhanden sein"""
        for n, rule_2_mere in enumerate(start_library[str(segment_number)]):  # 0. Spacer = 2. Segment
            if rule_2_mere == sequence[0: 2]:  # geht die Sachen einzeln durch und vergleicht diese. Wenn identisch, dann..
                decision_support += "T"  # True gemerkt
                break

        if len(decision_support) != 1:  # wenn noch nicht in die Liste hinzugefügt --> False merken
            decision_support += "F"

        """Restliche 2-mere nach 2. Tabelle überprüfen"""
        for o in range(0, len(sequence) - 2, 2):  # restliche Sequenz mit je 2 Basen durchgehen
            current_2_mere = sequence[o: o + 2]  # 1. 2-mere
            next_2_mere = sequence[o + 2: o + 4]  # nächste 2-mere
            counter_2mere += 1  # nächster 2-mere
            for j, rule_2_mere in enumerate(map_library[current_2_mere]):
                if rule_2_mere == next_2_mere:
                    decision_support += "T"
                    break

            if len(decision_support) != counter_2mere:  # Wenn noch nichts zur Liste hinzugefügt
                decision_support += "F"  # --> False merken

    else:
        counter_2mere = 0  # Zähler welcher 2-mere gerade betrachtet wird
        """Restliche 2-mere nach 2. Tabelle überprüfen"""
        for o in range(0, len(sequence) - 2, 2):  # restliche Sequenz mit je 2 Basen durchgehen
            current_2_mere = sequence[o: o + 2]  # 1. 2-mere
            next_2_mere = sequence[o + 2: o + 4]  # nächste 2-mere
            counter_2mere += 1  # nächster 2-mere
            for j, rule_2_mere in enumerate(map_library[current_2_mere]):
                if rule_2_mere == next_2_mere:
                    decision_support += "T"
                    break

            if len(decision_support) != counter_2mere:  # Wenn noch nichts zur Liste hinzugefügt
                decision_support += "F"  # --> False merken

    return decision_support


def translate_dna_to_rs_singular(segment, index):   # Will er nicht ausführen? Wieso RS? Wo ist der im Code?
    """
    DNA --> Binär --> AscII Latin 1, wo ist RS drin?
    :param segment:
    :param index:
    :return:
    """

    binary_segment = translate_dna_to_binary_singular(segment, index)  # DNA --> Binär mit DMR Schema

    # Binärzahlen in eine Liste mit Segmenten aus je 8 bits einteilen: 0b + 8 bits
    binary_string_8_bit_list = ["0b" + binary_segment[8 * i: 8 * (i + 1)] for i in range(len(binary_segment) // 8)]

    # Erzeugen einer Liste mit den Werten der gespeicherten Bits Segmente: Bsp.: int('10101010',2) = 170, (,2) = Base 2
    integer_value_list = [int(bit_string, 2) for bit_string in binary_string_8_bit_list]

    # Übersetzen der Werte in universelle Zeichencodierung Bsp: chr(97) = A und Liste zusammenfügen standard Ascii Liste
    ascii_rs_string = "".join([chr(integer_value) for integer_value in integer_value_list])

    # Ascii Latin 1 (andere Ascii Liste)
    ascii_latin1_rs_string = ascii_rs_string.encode("latin1")  # war einfacher erneut zu übersetzen

    return ascii_latin1_rs_string


def translate_decimal_to_dmr_singular(segment, index):
    """
    This function is to translate a segment with decimal numbers into DNA according to the DMR scheme.
    For this, the segment is first translated into binary numbers and then the DMR scheme is applied.
    :param segment: Input a segment with decimal numbers.
    :param index: Describes in which segment the DNA is and gives the DMR scheme the first coding.
    :return: Returns a DNA strand which has been translated using the DMR method.
    """
    segment_count = index % 4  # Beschreibt in welchem Segment DNA ist und gibt dem DMR Schema die erste Kodierung vor.
    binary_string = ""

    for value in segment:
        binary_string += format(value, '#010b')[2:]  # Wert für Wert durchgehen und übersetzen.
        # [2:]: print 0b + sequenz, 0b wollen wir nicht haben: also nur die Sequenz
    dna_string = translate_binary_to_dmr_singular(binary_string, segment_count)

    return dna_string



#############################
##### Check DMR Scheme ######
#############################


def determine_correct_segments(sequence_list, rs_codec_value, segment_count):
    """
    This Method takes a list of DNA sequences and checks if they are DMR and RS decodeable. To do that it translates
    the sequences to a decimal Format to RS decode it. Finally, the function returns all DMR sequences.

    Args:
        sequence_list: mutated(Substitution) DNA string
        rs_codec_value: The value of the RSCodec. In other words the number of decimal Values that have been added by
            Reed Solomon.
        segment_count:

    Returns:
         decoded_sequences_list: all possible options of the original DNA strand
    """
    rs_coder = RSCodec(rs_codec_value)

    decoded_sequences_list = []
    for i, DNA_error_segment in enumerate(sequence_list):
        validation_list = validate_list(DNA_error_segment, segment_count)
        validation_list_check = check_validation(validation_list)   # Ausgabe True oder False, wenn falsche vorhanden sind
        if validation_list_check:
            rs_segment = translate_dna_to_rs_singular(DNA_error_segment, segment_count)

            try:
                _, rs_string_decoded, errata_pos = rs_coder.decode(rs_segment)

                if rs_coder.check(rs_string_decoded):
                    dna_decoded_segment = translate_decimal_to_dmr_singular(rs_string_decoded, segment_count)
                    decoded_sequences_list += [(DNA_error_segment, dna_decoded_segment, i)]

            except:
                pass

    return decoded_sequences_list


def validate_list(segment, segment_count):
    """
    This function returns an array that has a corresponding label to each status the two-mere on the same index has
    in regard to the encoding rule.
    :param segment: Input of the segment of a DNA sequence.
    :param segment_count: Definition of the segment number of the entered DNA strand.
    :return: The function outputs a validation list. This shows whether a 2-mere matches the previous 2-mere in the
     DMR scheme. It also determines whether the initial 2-mere exists in the scheme.
    """

    # label:
    # correct does not mean that the two-mere is the same as the corresponding two-mere in the original string
    # it means that it fits into the DMR scheme

    """1. Abrufen der Mapping Libraries und Vorbereitung"""
    # Lade die JSON-Datei in ein Dictionary
    root_dir = os.path.abspath(os.curdir)

    with open(f"{root_dir}\\mapping_table_dmr.json") as f:
        data = json.load(f)
    map_library = data["map_library"]
    map_library_start_2_mere = data["initial_2mer"]


    # Eingegebene Sequenz wird in je 2 Basen segmentiert
    segmentated_string = segment_string(segment, 2)
    # Definieren erstes und letztes 2-mere
    start_2_mere = segmentated_string[0]
    last_2_mere = segmentated_string[-1]
    validation_list = []

    segment_mod = segment_count % 4

    """2. Validierung"""

    # <editor-fold desc="2.1 Validierung des Start 2-meres.">

    # start-two-mere
    # sT_nmT: start-two-mere --> correct | next-two-mere --> correct
    # sT_nmF: start-two-mere --> correct | next-two-mere --> false
    # sF: start-two-mere --> false

    if start_2_mere in map_library_start_2_mere[str(segment_mod)]:   # Falls das Start 2-mere in Map. Tab. 1 vorhanden

        if segmentated_string[1] in map_library[start_2_mere]:    # [1] das 2. 2-mere in 2. Mapping Tabelle von 1. mere
            validation_list += ["sT_nmT"]                         # beides da
        else:
            validation_list += ["sT_nmF"]       # erste richtig, zweite falsch
    else:
        validation_list += ["sF"]               # Start falsch
    # </editor-fold>

    # <editor-fold desc="2.2 Validierung der restlichen 2-mere.">

    # two-mere
    # tmT_nmT: two-mere --> correct | next-two-mere --> correct
    # tmT_nmF: two-mere --> correct | next-two-mere --> false
    # tmF_nmT: two-mere --> false | next-two-mere --> true
    # tmF_nmF: two-mere --> false | next-two-mere --> false

    # Bsp.:   AT  TC    GC:    i muss in der Tabelle von i-1 und i+1 in der Tabelle von i sein
    #        i-1   i   i+1
    for i in range(1, len(segmentated_string) - 1):                              # ohne erste und letzte 2-mere

        if segmentated_string[i] in map_library[segmentated_string[i - 1]]:      # i in Mapping Tabelle i+1

            if segmentated_string[i + 1] in map_library[segmentated_string[i]]:  # i+1 in Mapping Tabelle i
                validation_list += ["tmT_nmT"]

            else:
                validation_list += ["tmT_nmF"]                                   # i+1 nicht in Mapping Tabelle i

        else:                                                                    # i nicht in Mapping Tabelle i+1

            if segmentated_string[i + 1] in map_library[segmentated_string[i]]:  # i+1 in Mapping Tabelle i
                validation_list += ["tmF_nmT"]

            else:                                                                # i+1 nicht in Mapping Tabelle i
                validation_list += ["tmF_nmF"]

    # </editor-fold>

    # <editor-fold desc="2.3 Validierung des letzten 2-meres.">

    # last-mere
    # lT: last-two-mere --> correct
    # lF: last-two-mere --> false

    if last_2_mere in map_library[segmentated_string[-2]]:  # letzter in Mapping Tabelle vorletzter 2-mere
        validation_list += ["lT"]

    else:
        validation_list += ["lF"]
    # </editor-fold>

    return validation_list


def check_validation(validation_list):
    # sT_nmT: start-two-mere --> correct | next-two-mere --> correct
    # sT_nmF: start-two-mere --> correct | next-two-mere --> false
    # sF: start-two-mere --> false

    # two-mere
    # tmT_nmT: two-mere --> correct | next-two-mere --> correct
    # tmT_nmF: two-mere --> correct | next-two-mere --> false
    # tmF_nmT: two-mere --> false | next-two-mere --> true
    # tmF_nmF: two-mere --> false | next-two-mere --> false

    # last-mere
    # lT: last-two-mere --> correct
    # lF: last-two-mere --> false
    for entry in validation_list:
        if entry == "sT_nmF" or entry == "sF" or entry == "tmT_nmF" or entry == "tmF_nmT" or \
                entry == "tmF_nmF" or entry == "lF":
            return False

    return True


def find_inconsistency_indices_new(validation_list):
    """
    This function should detect all discrepancies in the validation list and output the positions which do not match
    the DMR scheme.
    :param validation_list: This list contains evaluations whether the first 2-mere, the following 2-mere and the last
    2-mere fit into the DMR scheme. For this purpose, the current 2-mere is compared with the one before it and
    additionally the 2-mere following it is compared with the current one and evaluated whether both variants fit
    into the scheme.
    :return: The function returns an array containing locations that do not fit the DMR scheme.
    """

    inconsistency_index_array = []
    # Wo befinden sich diese Abkürzungen, welche für fehlerhafte Stellen stehen
    for i, status in enumerate(validation_list):
        if status == "sF" or status == "sT_nmF" or status == "tmT_nmF" or status == "tmF_nmF" or status == "tmF_nmT" or status == "lF":
            inconsistency_index_array += [i]

    return inconsistency_index_array


def get_inconsistency_neighbours_new(validation_list):
    """
       Determines the size of consecutive Inconsistencies and returns them as a list.

       Args:
           validation_list: Validation list of DMR Sequence

       Returns:
            consecutive_inconsistencies_list: list that contains the index and size of an inconsistency chain
    """
    inconsistency_neighbour_list = []
    inconsistency_neighbours = []
    switch = False

    for i, status in enumerate(validation_list):
        if status == "sF" or status == "sT_nmF" or status == "tmT_nmF" or status == "tmF_nmF" or status == "tmF_nmT" or status == "lF":

            if switch == False:

                switch = True
                inconsistency_neighbours += [i]

            else:
                inconsistency_neighbours += [i]

        else:
            if switch == True and len(inconsistency_neighbours) > 1:
                inconsistency_neighbour_list += [inconsistency_neighbours]
                inconsistency_neighbours = []
                switch = False

            else:
                inconsistency_neighbours = []
                switch = False

    if len(inconsistency_neighbours) != 0:
        inconsistency_neighbour_list += [inconsistency_neighbours]
        inconsistency_neighbours = []

    return inconsistency_neighbour_list


def check_for_neighboring_inconsistency(validation_list, position):
    """
    This function compares two places of the validation list with each other.
    If both entries do not match the DMR schema, then True is returned, otherwise false.
    :param validation_list: This list contains evaluations whether the first 2-mere, the following 2-mere and the last
    2-mere fit into the DMR scheme. For this purpose, the current 2-mere is compared with the one before it and
    additionally the 2-mere following it is compared with the current one and evaluated whether both variants fit
    into the scheme.
    :param position: Defines a position of the validation list which should be checked.
    :return: True --> Both positions are erroneous. False --> Only one position is erroneous.
    """

    # <editor-fold desc="1. Anschauen des ersten Eintrages der Validation Liste">

    if position == 0:
        status_1 = validation_list[position]        # erster Eintrag
        status_2 = validation_list[position + 1]    # zweite Eintrag

        if status_1 == "lF" or status_1 == "tmF_nmF" or status_1 == "tmF_nmT" or status_1 == "sF":

            if status_2 == "lF" or status_2 == "tmF_nmF" or status_2 == "tmF_nmT" or status_2 == "sF":

                return True     # 1. und 2. Einträge enthalten fehlerhafte Einträge

        else:
            return False        # nur der 1. Eintrag enthält einen fehlerhaften Eintrag
    # </editor-fold>

    # <editor-fold desc="2. Anschauen des letzten Eintrages der Validation Liste">

    elif position == len(validation_list) - 1:
        status_1 = validation_list[position]        # letzter Eintrag
        status_2 = validation_list[position - 1]    # vorletzter Eintrag

        if status_1 == "lF" or status_1 == "tmF_nmF" or status_1 == "tmF_nmT" or status_1 == "sF":

            if status_2 == "lF" or status_2 == "tmF_nmF" or status_2 == "tmF_nmT" or status_2 == "sF":
                return True

        else:
            return False
    # </editor-fold>

    # <editor-fold desc="3. Anschauen einer der restlichen Einträge der Validation Liste">

    else:
        status_1 = validation_list[position]
        status_2 = validation_list[position - 1]
        status_3 = validation_list[position + 1]

        if status_1 == "lF" or status_1 == "tmF_nmF" or status_1 == "tmF_nmT" or status_1 == "sF":

            if status_2 == "lF" or status_2 == "tmF_nmF" or status_2 == "tmF_nmT" or status_2 == "sF":
                return True                # Position + Postion -1 fehlerhaft

            elif status_3 == "lF" or status_3 == "tmF_nmF" or status_3 == "tmF_nmT" or status_3 == "sF":
                return True                # Position + Postion +1 fehlerhaft

        else:
            return False
    # </editor-fold>



######################
##### Encode RS ######
######################


def rs_encode_with_symbol_number(data: list, symbol_number: int, exp: int):
    """
    Adds RS codes to a set of bytes objects at a defined maximum percent error able to be corrected

    Args:
        data: list in which the values are bytes objects to which RS codes are added. The List can also contain only
              one item but has to be a list in a list. [[1,2,3]]
        symbol_number: Number of RS Symbols to be added.
        exp: exponent which defines the size of the Galois field. For example, exp = 12 defines a GF of 4095 (2^12 - 1)

    Returns:
        encoded: List with the bytes objects and the added RS Codes.
        max_errors: Returns the number of errors that can be corrected.
        max_erasures: Returns the number of erasures that can be corrected.
    """
    # Determine which RS Codec to use for each based on the percent error

    encoded = []
    rsc = RSCodec(symbol_number, c_exp=exp)  # adding a defined number of RS symbols

    for segment in data:
        if len(data) == 1:
            encoded = rsc.encode(segment)
        else:
            encoded.append(list(rsc.encode(segment)))

    max_errors, max_erasures = rsc.maxerrata(verbose=False)

    return encoded, max_errors, max_erasures


def bits_to_dna(fp_root="", *args: str):
    """
    Translates strings of bits (0 or 1) into strings of DNA bases (A, C, T, or G) using an alternative method

    Args:
        fp_root: Possible addition of filepath needed for demonstrator. For other script just use it by setting this "".
        *args: strings of bit sequences that need to be translated to DNA

    Returns:
        Translated DNA strings in the order that the corresponding bit strings were inputted.
    """

    with open(fp_root + "mapping_table_two_bit.json") as json_file:
        mapping = json.load(json_file)["map_library"]

    dna_sequences = []
    for bit_string in args:
        two_bits_list = [bit_string[i:i+2] for i in range(0, len(bit_string), 2)]
        dna_sequences.append(''.join([random.choice(mapping.get(tb)) for tb in two_bits_list]))


    return tuple(dna_sequences)


######################
##### Decode RS ######
######################


def rs_force_uncorrect(codec_size: int, exp: int, *args):
    """
    Removes the RS codes from a list of bytes without correcting the string
    Args:
        codec_size: size of the codec
        exp: exponent which defines the size of the Galois field. For example, exp = 12 defines a GF of 4095 (2^12 - 1)
        **args: list of byte lists to decode
    Returns:
        The decoded variables in the order in which they were called in the function
    """

    decoded = []
    for bitlist in args:
        segment_length = 2**exp-1
        splits = [bitlist[i:i+segment_length] for i in range(0, len(bitlist), segment_length)]
        uncorrected_chunked = list([chunk[:len(chunk)-codec_size] for chunk in splits])
        uncorrected_flat = [item for sublist in uncorrected_chunked for item in sublist]
        decoded.append(uncorrected_flat)

    return tuple(decoded)


def dna_to_bits(*args: str):
    """
    Translates strings of DNA (A, C, T, or G) into strings of bits (0 or 1) using a mapping table.

    Args:
        *args: strings of DNA sequences that need to be translated to bits

    Returns:
        Translated bit strings in the order that the corresponding DNA strings were inputted.
    """


    with open("mapping_table_two_bit.json") as json_file:
        old_mapping = json.load(json_file)["map_library"]

    vals = list(old_mapping.values())
    keys = list(old_mapping.keys())

    mapping = {}
    for i in range(len(vals)):
        for val in vals[i]:
            mapping[val] = keys[i]

    bit_sequences = []
    for bit_string in args:

        two_bits_list = [bit_string[i:i+2] for i in range(0, len(bit_string), 2)]
        bit_sequences.append(''.join([mapping.get(tb) for tb in two_bits_list]))

    return tuple(bit_sequences)


##########################
##### Random Errors ######
##########################


def binom_mutations_with_spacer_ignorance(sequence: str, dna_length_without_X: int, substitution_freq: float, insertion_freq: float, deletion_freq: float,
                                          binom=True):
    """
    A function used to generate random substitutions, insertions, and deletions in a DNA string given the DNA string and the frequency of each event.
    The DNA can contain spacers, which are marked with XXXXX. However, the errors are inserted only in the DNA sequence and not in the spacer sequences,
    so that these can be removed later with split more easily.Also, able to specify whether the event frequencies are the probability of a single
    event (meaning the total number of events follows a binomial distribution), or whether the event frequency is a percent of the length of the sequence.

    Args:
        sequence: DNA sequence string
        dna_length_without_X: Input of the DNA length without the spacer sequences. This is needed to calculate the number of individual errors.
        substitution_freq: decimal between 0 and 1.0
        insertion_freq: decimal between 0 and 1.0
        deletion_freq: decimal between 0 and 1.0
        binom: boolean value specifying whether the total number of events for each event is calculated following a
               binomial distribution or not

    Returns:
        DNA or binary sequence string with the number of events given by the event frequencies.
        The locations of the events are determined by a random sample of the list of bases or bits.
    """
    # check language
    if set(np.unique(list(sequence))).issubset({'G', 'C', 'A', 'T', 'X'}):
        alphabet = ["A", "C", "T", "G"]
    elif set(np.unique(list(sequence))).issubset({'0', '1'}):
        alphabet = ["0", "1"]
    else:
        raise ValueError(
            "Incorrect sequence language. Make sure your sequence only includes '0' or '1' if it is a binary sequence, "
            "or 'G', 'C', 'A', 'T' or 'X' if it is a DNA sequence")

    # set number of events
    if binom:
        num_to_substitute = stats.binom.rvs(dna_length_without_X, substitution_freq)
        num_to_insert = stats.binom.rvs(dna_length_without_X, insertion_freq)
        num_to_delete = stats.binom.rvs(dna_length_without_X, deletion_freq)
    else:
        num_to_substitute = int(dna_length_without_X * substitution_freq)
        num_to_insert = int(dna_length_without_X * insertion_freq)
        num_to_delete = int(dna_length_without_X * deletion_freq)

    # substitutions
    # substitute = random.sample(range(len(sequence)), k=num_to_substitute)
    substitute = random.sample([pos for pos, letter in enumerate(sequence) if letter in alphabet], k=num_to_substitute)
    sequence = [
        (random.choice([base for base in alphabet if base != (sequence[i])]) if i in substitute else sequence[i]) for i
        in list(range(len(sequence)))]
    sequence = ''.join(sequence)

    # insertions
    insert = random.sample([pos for pos, letter in enumerate(sequence) if letter in alphabet], k=num_to_insert)
    sequence = [(sequence[i] + random.choice(alphabet) if i in insert else sequence[i]) for i in
                list(range(len(sequence)))]
    sequence = ''.join(sequence)

    # deletions
    delete = random.sample([pos for pos, letter in enumerate(sequence) if letter in alphabet], k=num_to_delete)
    sequence = [sequence[i] for i in list(range(len(sequence))) if i not in delete]
    sequence = ''.join(sequence)

    return sequence



def single_substitution_correction(segment, error_two_mere_position):
    """
    This Method returns all possible options of an original DNA strand before it had a single substitution error
    introduced at a certain position.

    Args:
        error_two_mere_position:
        segment:

    Returns:
        correction_option_list: all possible options of the original DNA strand
    """

    first_base_position = error_two_mere_position * 2
    all_two_meres = ["AA", "CC", "GG", "TT", "AC", "CG", "GT", "TA", "AG", "CT", "GA", "TC", "AT", "CA", "GC", "TG"]

    segment_options = []

    for two_mere in all_two_meres:
        new_segment = segment[:first_base_position] + two_mere + segment[first_base_position + 2:]
        segment_options += [new_segment]

    return segment_options, [error_two_mere_position * 2, "M"]

#########################
##### Edit Distance #####
#########################


def score_pairs_fast(pairs: list, as_percent=True):
    """
    A function that calculates the Levenshtein distance between two strings, used here for binary or DNA sequences

    Args:
        as_percent:
        pairs: A list of tuples, where each tuple has 3 strings. The first two elements that are compared to each other,
               and the second is the name for the pair
               ie [(string1, string2, name1_2), (string3, string4, name3_4)]

    Returns:
        A list of tuples where each tuple has the format (name1_2, score1_2). The score is calculated by first computing
        the Levenshtein distance between the two strings and then dividing by the length of the longest string. This
        value is then subtracted from 1 to yield a percent similarity from [0, 1]. Uses teh Levenshtein package
        ex: score_pairs([("kitten", "sitting", "kitten and sitting"), ("Sunday", "Saturday", "Sunday and Saturday")])
        returns  [('kitten and sitting', 0.5714285714285714), ('Sunday and Saturday', 0.625)]
    """
    similarity_scores = []
    names = []

    for pair in pairs:
        s1 = str(pair[0])
        s2 = str(pair[1])
        name = str(pair[2])

        if len(s1) == 0 or len(s2) == 0:
            edit_dist = max(len(s1), len(s2))
        else:
            edit_dist = Levenshtein.distance(s1, s2)

        if as_percent:
            similarity_scores.append(1 - edit_dist / max(len(s1), len(s2)))
        else:
            similarity_scores.append(edit_dist)
        names.append(name)

    return list(zip(names, similarity_scores))
