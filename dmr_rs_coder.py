import os
import json
import math
import warnings
import functions as func
from reedsolo import RSCodec
import dmr_level_master


class DMR_RS_Coder:
    """
    This Class is a decoder/encoder in similar style to the Reed-Solomon decoder/encoder.
    The class (in construction) should contain Methods to encode Binary to DNA and decode DNA to Binary.
    The used Mapping-Rule DMR (Dynamic Mapping Rule) is an attempt to utilize a self-correcting DNA code
    Inbetween a Reed-Solomon encoding is used to further enhance the correctability of the sequence. Reed-Solomon and
    DMR interact with each other so that (in construction) in synergy a better correction should be possible.

    Args:
       rs_codec_value: The values of the RSCodec of the used Reed-Solomon code
       segment_size: The size of the DMR segment
       mapping_rule: The DMR encoding rule
       --> further described in dmr_translator.py
       mapping_rule_start_2_mere : The DMR encoding rule that is utilized for the first two encoded Bases of a segment
       --> further described in dmr_translator.py

       min_codec_value : The minimum RSCodec value per segment. Since the codec is downsized to its absolute minimum
                         this Parameter can set a boundary so that the Codec value isn't downsized beyond this point
                         --> see recalculate_coded ()

    Returns:
          xxxxxxxxxxxxxxxxxxxxxx
    """

    def __init__(self, rs_codec_value: int, min_codec_value: int = 0, min_segment_length: int = 0):

        # Lade die JSON-Datei in ein Dictionary
        root_dir = os.path.abspath(os.curdir)
        with open(f"{root_dir}\\mapping_table_dmr.json") as f:
            data = json.load(f)
        # Zugriff auf den jeweiligen Teil des Dictionaries
        map_library = data["map_library"]
        map_library_start_2_mere = data["initial_2mer"]

        self.rs_codec_value = rs_codec_value  # Beschreibt den normalen Codec Value bevor er auf die verkleinerten Segment angepasst wird (Bsp. RS-CODEC 32)
        self.min_codec_value = min_codec_value  # Beschreibt den eingestellen runtergerechneten Codec Value der die verkleinerten Segmente angepasst wurde
        self.min_segment_length = min_segment_length  # Beschreibt die eingestellte runtergerechnete Segmentlänge
        self.mapping_rule = map_library
        self.mapping_rule_start_2_mere = map_library_start_2_mere
        self.RS_new_codec_size = None

        if rs_codec_value < min_codec_value:
            warnings.warn(f"The minimal downsized Codec Value is bigger than the defined codec value!!! The minimal "
                          f"codec value is denied and wont be taken into account")
            self.min_codec_value = 0


    def recalculate_codec(self):
        """
            To correct the DMR the Reed-Solomon Codec should be as small as possible
            The segment size is of DMR is as Big as the Reed-Solomon segments in bits (RS-segment-size * 8)
            It´s harder to correct many errors in one big segment than a few errors in many small segments
            Therefore the smallest possible Reed-Solomon size at which no error correction utility is lost is calculated
            It is important to be noted that the calculated codec size is fixed.
            This is the used RSCodec sized

            Returns:
                    s_new_codec_size, payload_size
        """

        payload_size = 0
        rs_new_codec_size = 1

        # Hier wird der Codec/Segmentlänge entsprechend der angegebenen Parameter runtergerechnet ( min_codec_value, min_segment_length)
        # Wenn die Parameter sich wiedersprechen wird der minimale Codec an die minimale Segmentlänge angepasst

        if self.min_segment_length == 0:

            if self.min_codec_value == 0:
                payload_size = math.floor((255 - self.rs_codec_value) / self.rs_codec_value)

            else:

                payload_size = math.floor(self.min_codec_value * (255 - self.rs_codec_value) / self.rs_codec_value)

                if payload_size != 0:
                    rs_new_codec_size = self.min_codec_value

                else:

                    maximal_downsized_codec = 1 / ((255 - self.rs_codec_value) / self.rs_codec_value)
                    warnings.warn(f"The minimal downsized Codec Value is too big for the defined codec value (payload "
                                  f"would be zero in length)!!! The minimal codec value has been lowered"
                                  f" {self.min_codec_value} --> {maximal_downsized_codec}")
                    payload_size = math.floor(maximal_downsized_codec * (255 - self.rs_codec_value) /
                                              self.rs_codec_value)

        else:
            if self.min_codec_value == 0:
                rs_new_codec_size = math.ceil(self.min_segment_length * (self.rs_codec_value / 255))
                payload_size = self.min_segment_length - rs_new_codec_size

                if rs_new_codec_size == 1 and self.min_segment_length == 1:
                    payload_size = 1

            else:
                maximal_downsized_codec = math.ceil(self.min_segment_length * (self.rs_codec_value / 255))

                if self.min_codec_value < maximal_downsized_codec:
                    rs_new_codec_size = maximal_downsized_codec
                    payload_size = self.min_segment_length - rs_new_codec_size

                    warnings.warn(f"The minimal Codec Value is too small for the defined minimal segment length!!! The "
                                  f"minimal codec value has been raised {self.min_codec_value} --> {rs_new_codec_size}")

                else:
                    rs_new_codec_size = self.min_codec_value
                    payload_size = self.min_segment_length - maximal_downsized_codec

        self.RS_new_codec_size = rs_new_codec_size
        return rs_new_codec_size, payload_size


    def encode_rsm(self, data):
        """
        To correct the DMR the Reed-Solomon Codec should be as small as possible
        The segment size is of DMR is as Big as the Reed-Solomon segments in bits (RS-segment-size * 8)
        It´s harder to correct many errors in one big segment than a few errors in many small segments
        Therefore the smallest possible Reed-Solomon size at which no error correction utility is lost is calculated
        It is important to be noted that the calculated codec size is fixed.
        This is the used RSCodec sized

        Args:
                data: The data that is encoded to Reed-Solomon code

        Returns:
                 rsm_segments, new_codec_size, payload_size
        """
        new_codec_size, payload_size = self.recalculate_codec()  # Hier wird die neue Codec Größe in der Funktion recalculate_codec eingestellt und zusätzlich
        # werden der neue Codec und Segmentlänge zurückgegeben
        rs_coder = RSCodec(new_codec_size)

        payload_segments = func.segment_string(data, payload_size)  # Segmentation der DNA

        rsm_segments = [rs_coder.encode(segment) for segment in payload_segments]  # Reed-Solomon Codierung

        return rsm_segments, new_codec_size, payload_size


    def initial_scan_correction(self, sequence_list, index=None):
        """
        This Function attempts correction of the given DMR Segmented Fragments by Reed-Solomon only. This serves as a first scan to minimize DNA segments that
        need to be corrected by DMR

        Args:
                sequence_list: Segmented DNA list of DMR-RSM encoded segments
                index: If given all entries of the sequence list are attempted to be corrected with the given index DMR Start Two-mere rule

        Returns:
                 erroneous_segment_list, decoded_segment_list
        """
        rs_coder = RSCodec(self.RS_new_codec_size)

        decoded_segment_list = []

        erroneous_segment_list = []

        for i, segment in enumerate(sequence_list):
            if index == None:
                index_choice = i
            else:
                index_choice = index

            validation_list = func.validate_list(segment, index_choice)
            validation_list_check = func.check_validation(validation_list)
            # Validation List is initialized to check for inconsistencies in the DMR scheme

            if validation_list_check:
                rs_segment = func.translate_dna_to_rs_singular(segment, index_choice)
                # If the segment passes the validation list test (no inconsistencies in segment) it is attempted to be RS decoded
                try:
                    _, rs_string_decoded, errata_pos = rs_coder.decode(rs_segment)

                    if rs_coder.check(rs_string_decoded):
                        decoded_segment_list += [(segment, index_choice)]

                    else:
                        erroneous_segment_list += [(segment, index_choice)]

                except:
                    erroneous_segment_list += [(segment, index_choice)]

            else:
                erroneous_segment_list += [(segment, index_choice)]
            # All segments that were decodeable get added to a decodeable segment list
            # All segments that were not decodeable get added to an erroneous segment list
        return erroneous_segment_list, decoded_segment_list


    def after_scan_correction_try_Jess(self, correction_list, verbose=True):
        """
        This Function attempts correction of the given DMR Segmented Fragments with the DMR Scheme.
        Args:
            correction_list: Segmented DNA list of DMR-RSM encoded segments
            verbose: Determine if all results should be displayed.
        Returns:
            errornous_sequences, corrected_sequences
        """
        corrected_sequences = []
        errornous_sequences = []

        for sequence_sub_list in correction_list:
            error_sequence, sequence_count = sequence_sub_list
            sequence_list, level, segment_count = dmr_level_master.soft_force_new(error_sequence, sequence_count, self.RS_new_codec_size, verbose=verbose)
            # Diese Funktion nimmt eine Sequenz und die bekannten Parameter (Segment_Index, Codec) und versucht diese durch DMR zu korrigieren

            # Liste mit Möglichkeiten aller korrigierbarer Sequenzen --> oft gleich --> wenn nicht gleich dann oft false positiv.
            # Er geht alle Level bis etwas korrigierbares gefunden wurde durch und schaut dann die Ergebnisse an
            corrected_sequence_option_list = []
            for liste in sequence_list:
                corrected_sequence_option_list += [(liste[1], segment_count)]

            switch = True
            # Switch: Wenn alle Sequenzen gleich sind, dann war die Sequenz erfolgreich dekodierbar.
            # Wenn eine Sequenz unterschiedlich ist, dann kommt das in die nicht dekodierbare Liste.
            for entry in corrected_sequence_option_list:
                if corrected_sequence_option_list[0][0] != entry[0]:
                    switch = False

            # Ausprobieren: versuchen aus den gefundenen Sequenzen die richtige zu finden!
            if switch == False:

                individual_sequences = list(set(corrected_sequence_option_list))
                individual_sequences_count = [corrected_sequence_option_list.count(entry) for entry in individual_sequences]
                max_seq = max(individual_sequences_count)
                max_count = individual_sequences_count.count(max_seq)

                # Mehrheitsentscheid gleicher Sequenzen
                if max_count == 1:
                    corrected_sequence_option_list = [individual_sequences[individual_sequences_count.index(max_seq)]]
                    switch = True

                # Auswahl nach Editdistanz
                else:
                    edit_distance = [(func.score_pairs_fast([(error_sequence, item[0], "DNA")], as_percent=True), item) for item in individual_sequences]
                    max_edit, max_edit_sequence = max(edit_distance)
                    corrected_sequence_option_list = [max_edit_sequence]
                    switch = True

                # Mehrheitsentscheid gleicher Sequenzabschnitte --> Erstellung einer Sequenz daraus und überprüfen, ob in Liste vorhanden.

            if switch == False or len(sequence_list) == 0:
                errornous_sequences += [(error_sequence, sequence_count)]
            elif switch and len(corrected_sequence_option_list) > 0:
                corrected_sequences += [corrected_sequence_option_list[0]]
            else:
                corrected_sequences += corrected_sequence_option_list[0]

        return errornous_sequences, corrected_sequences
