import os
import json
import functions as func


def soft_force_new(error_segment, segment_count, rs_codec_value, verbose=True):
    """
    This Method tries to correct DMR Segments in a rather low level approach. Depending on the difference of original
    segment length and known segment length (probably only one I don´t think I'll do more than that) this Method is used
    It only tries to correct with a basic approach without expecting complicated error patterns like multiple Insertions
    /Deletions that leave the message errourness but at its original size.
    The errors this Method tries to correct are:
    -single Deletion Errors with and without single Substitution Error
    -single Insertion Errors with and without single Substitution Error
    -Substitution Errors (I still have to decide whether we will do only 1 or 2 Substitution errors and also undetected ones --> I guess yes lets see)

    Args:
        rs_codec_value:
        error_segment: errourness DNA string
        segment_count: Number of point mutations to be introduced
        verbose: Determine if all results should be displayed.
    Returns:
        validation_list = DMR_translator.validate_list(segment,segment_count)
    """


    level = 0

    # für später ausprobieren --> Insertion oder Deletion vorhanden? Gehe gleich zu Level 5
    # if len(error_segment) % 8 != 0:
    # level = 5

    if level == 0:
        correction = soft_level_0(error_segment, segment_count, rs_codec_value)

        if correction:
            if verbose:
                print(f"Segmentcount: {segment_count}, Level: {level}")
            return correction, level, segment_count

        else:
            level += 1

    if level == 1:
        correction = soft_level_1(error_segment, segment_count, rs_codec_value)

        if correction:
            if verbose:
                print(f"Segmentcount: {segment_count}, Level: {level}")
            return correction, level, segment_count

        else:
            level += 1

    if level == 2:
        correction = soft_level_2(error_segment, segment_count, rs_codec_value)

        if correction:
            if verbose:
                print(f"Segmentcount: {segment_count}, Level: {level}")
            return correction, level, segment_count

        else:
            level += 1

    if level == 3:
        correction = soft_level_3(error_segment, segment_count, rs_codec_value)

        if correction:
            if verbose:
                print(f"Segmentcount: {segment_count}, Level: {level}")
            return correction, level, segment_count

        else:
            level += 1


    return correction, level, segment_count


""" 
Level_Description: Beim Nachvollziehen der Funktion die Levelbeschreibungen zusätzlich angucken

soft_level_0
Es wird davon ausgegangen, dass nur Substitutionen vorhanden sind. Sind Inkonsistenzen vorhanden so wird das Mapping Schema genutzt, um diese mit den
vorhandenen Möglichkeiten zu korrigieren, wobei mindestens 1 Base zum vorgegebenen 2-mere passen muss. Ein 2-mere mit 2 falschen Basen kann hier 
nur bei dem Fehlermuster 'tmT_nmF', 'tmF_nmT' korrigiert werden. 
--> Erstes und letztes 2-mere: Hier wird das Start-2-mere anhand der Starttabelle und der Segmentnummer bestimmt! 
    Das letzte 2-mere kann nicht durch das darauffolgende korrigiert werden.  
--> 2 mere mittendrin: Möglichkeiten des vorherigen 2-meres mit passenden 2-meren anschauen und mit dem darauffolgenden vergleichen. 
--> 2 2-mere hintereinander: Brücke bauen und versuchen passende einzugrenzen.

soft_level_1
Es wird davon ausgegangen, dass nur Substitutionen vorhanden sind. Sind Inkonsistenzen vorhanden so wird das Mapping Schema genutzt, um diese mit den
vorhandenen Möglichkeiten zu korrigieren. Die Korrektur beruht auf dem von Level 0, aber ohne die Einschränkung der notwendigen passenden Base pro 2-mere. 
Zusätzlich wurde die Korrektur von 5 hintereinander folgenden falschen 2-meren dahingehend geändert, dass die 3 mittleren alle als falsch angenommen werden. 

soft_level_2
It is assumed that only substitutions are present and that no correction with level 0 and 1 has worked. Here, the locations of the inconsistencies are first 
determined without separating them as in levels 0 and 1. Due to substitutions, errors could occur that fit into the mapping scheme and thus output false 
inconsistencies. In this case, a correction should be made by attempting to correct the list with the longest inconsistency. This involves going through the 
2-mers and making individual substitutions by substituting each possible 2-mere (16 variants in total) at the position currently under consideration. 
The possibilities are corrected in level 0. If decodeable segments can be found, these are returned. If not, it continues with the next position. If no 
correctable segments are found with the substitution and level 0, the whole process is repeated with substitutions and a correction in level 1. 

"""


def soft_level_0(error_segment, segment_count, rs_codec_value):
    """
    It is assumed that only substitutions are present. If inconsistencies are present, the mapping scheme is used to correct these with the available options,
    whereby at least 1 base must match the specified 2-mere. A 2-mere with 2 incorrect bases can only be corrected here with the error pattern 'tmT_nmF',
    'tmF_nmT'. Level should correct substitutions that are recognised by the inconsistent indices.
    Error patterns are separated so that the positions can be considered individually!
    --> Separation of the error patterns according to 'tmF_nmT'.

    Args:
        error_segment: Input the errornous segment to correct.
        segment_count: Input the segment number of the segment for the DMR translation scheme.
        rs_codec_value: Input the used codec value for determination of the correct segments.

    Returns:
        List with the decoded sequence with were checked using the check_dmr_scheme.determine_correct_segments(possible_sequences, rs_codec_value,
        segment_count) function.
    """

    validation_list = func.validate_list(error_segment, segment_count)
    inconsistency_indices = func.find_inconsistency_indices_new(validation_list)
    inconsistency_count = len(inconsistency_indices)
    neighbouring_inconsistencies_list = func.get_inconsistency_neighbours_new(validation_list)
    len_seg = len(error_segment)
    num_mere = len_seg / 2

    for i, items in enumerate(neighbouring_inconsistencies_list):

        split_point = ['tmF_nmT']  # Definiere die Stellen, an denen die Liste geteilt werden soll
        result_list = []  # Initialisiere eine leere Liste für die geteilten Teillisten
        current_sublist = []  # Initialisiere eine Variable zum Zwischenspeichern der aktuellen Teilliste
        val = validation_list[items[0]:items[-1] + 1]
        for num, element in enumerate(val):
            # Füge das aktuelle Element zur aktuellen Teilliste hinzu
            current_sublist.append(items[num])

            # Überprüfe, ob das aktuelle Element eine der Split-Stellen ist
            if element in split_point:
                result_list.append(current_sublist)
                current_sublist = []

        if result_list != []:

            if current_sublist != []:
                result_list.append(current_sublist)

            neighbouring_inconsistencies_list.pop(i)
            for a in range(len(result_list), 0, -1):
                neighbouring_inconsistencies_list.insert(i, result_list[a - 1])

    len_neighbouring = [len(items) for items in neighbouring_inconsistencies_list]
    max_neighbour = max(len_neighbouring) if len_neighbouring else 0

    # keine Insertion/Deletion vorhanden, es müssen Fehler vom Mapping Schema da sein, höchstens 5 Inkonsistenzen hintereinander
    if inconsistency_count >= 1 and len(error_segment) % 2 == 0 and max_neighbour <= 5:

        #####################
        ### Vorbereitung ####
        #####################
        root_dir = os.path.abspath(os.curdir)
        with open(f"{root_dir}\\mapping_table_dmr.json") as f:
            data = json.load(f)
        map_library_start_2_mere = data["initial_2mer"]
        list_of_start_value = list(map_library_start_2_mere.values())
        map_library = data["map_library"]
        list_of_key = list(map_library.keys())
        list_of_value = list(map_library.values())
        trouble_parts = []

        for item in neighbouring_inconsistencies_list:
            last_2_mere = error_segment[item[0] * 2:item[0] * 2 + 2]
            current_2_mere = error_segment[item[1] * 2:item[1] * 2 + 2]
            possibilities = []

            ############################
            ### Start 2-mere falsch ####
            ############################
            # Start: F T
            if item == [0, 1] and validation_list[0: 2] == ['sF', 'tmF_nmT']:

                start_possibilities = list_of_start_value[segment_count % 4]
                for i, mere in enumerate(list_of_value):
                    if error_segment[2:4] in mere and error_segment[0] == list_of_key[i][0] and list_of_key[i] in start_possibilities or \
                            error_segment[2:4] in mere and error_segment[1] == list_of_key[i][1] and list_of_key[
                        i] in start_possibilities:  # Level 1 diese Zeile wegnehmen
                        possibilities.append(list_of_key[i] + error_segment[2:4])

                if possibilities != []:
                    trouble_parts.append([item, "First 2-mere: F, T", possibilities])


            # Start: F F T
            elif item == [0, 1, 2] and validation_list[0: 3] == ['sF', 'tmF_nmF', 'tmF_nmT']:
                next_2_mere = error_segment[item[2] * 2:item[2] * 2 + 2]

                # Variationen 1-mere
                variation_option_1 = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in list_of_start_value[segment_count % 4]:
                    if variation[0] == last_2_mere[0] or variation[1] == last_2_mere[1]:  # Level 1 diese Zeile wegnehmen
                        variation_option_1.append(variation)

                # Variationen 2-mere
                variation_option_2 = []
                for i, items in enumerate(list_of_value):
                    if next_2_mere in items:
                        if list_of_key[i][0] == current_2_mere[0] or list_of_key[i][1] == current_2_mere[1]:  # Level 1 diese Zeile wegnehmen
                            variation_option_2.append(list_of_key[i])

                # Variation 1 in Variation 2?
                possibilities = []
                for items in variation_option_1:
                    for j, rule_2_mere in enumerate(map_library[items]):
                        if rule_2_mere in variation_option_2:
                            possibilities.append(items + rule_2_mere + next_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "First 2-mere: F, F, T", possibilities])


            # 'sF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'
            elif item == [0, 1, 2, 3] and validation_list[0: 4] == ['sF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT']:

                errorneous_part_1 = error_segment[item[1] * 2:item[1] * 2 + 2]
                errorneous_part_2 = error_segment[item[2] * 2:item[2] * 2 + 2]
                third_2_mere = error_segment[item[3] * 2:item[3] * 2 + 2]

                # möglichen Start suchen
                start = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in list_of_start_value[segment_count % 4]:
                    if variation[0] == error_segment[0] or variation[1] == error_segment[1]:  # Level 1 diese Zeile wegnehmen
                        start.append(variation)

                # Variationen 1-mere
                variation_option_1 = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for items in start:
                    for variation in map_library[items]:
                        if variation[0] == errorneous_part_1[0] or variation[1] == errorneous_part_1[1]:  # Level 1 diese Zeile wegnehmen
                            variation_option_1.append(items + variation)

                # Variationen 2-mere
                variation_option_2 = []
                for i, items in enumerate(list_of_value):
                    if third_2_mere in items:
                        if list_of_key[i][0] == errorneous_part_2[0] or list_of_key[i][1] == errorneous_part_2[1]:  # Level 1 diese Zeile wegnehmen
                            variation_option_2.append(list_of_key[i])

                # Variation 1 in Variation 2?
                possibilities = []
                for items in variation_option_1:
                    for j, rule_2_mere in enumerate(map_library[items[2:]]):
                        if rule_2_mere in variation_option_2:
                            possibilities.append(items + rule_2_mere + third_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "First 2-mere: F, F, F, T", possibilities])


            ##############################
            ### letztes 2-mere falsch ####
            ##############################
            # 'tmT_nmF', 'tF'
            elif item == [num_mere - 2, num_mere - 1] and validation_list[-2:] == ['tmT_nmF', 'lF']:
                for variation in map_library[last_2_mere]:
                    if variation[0] == error_segment[-2] or variation[1] == error_segment[-1]:  # Level 1 diese Zeile wegnehmen
                        possibilities.append(last_2_mere + variation)

                if possibilities != []:
                    trouble_parts.append([item, "Last 2-mere: 'tmT_nmF', 'lF'", possibilities])


            # 'tmT_nmF', 'tmF_nmF', 'tF'
            elif item == [num_mere - 3, num_mere - 2, num_mere - 1] and validation_list[-3:] == ['tmT_nmF', 'tmF_nmF', 'lF']:
                next_2_mere = error_segment[item[2] * 2:item[2] * 2 + 2]
                # T F T
                for variation in map_library[last_2_mere]:
                    for j, rule_2_mere in enumerate(map_library[variation]):
                        # Rule in Variation stimmt überein → rule_2_mere als current_2_mere annehmen
                        if rule_2_mere == next_2_mere:  # variation ist der eigentliche 2-mere
                            # Variation aussortieren --> nur eine Substitution angenommen --> eine Base sollte passen!
                            if variation[0] == current_2_mere[0] or variation[1] == current_2_mere[1]:  # Level 1 diese Zeile wegnehmen
                                possibilities.append(last_2_mere + variation + next_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "Last 2-mere: 'tmT_nmF', 'tmF_nmF', 'lF'", possibilities])


            # 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'lF'  T, F, T, F
            elif item == [num_mere - 4, num_mere - 3, num_mere - 2, num_mere - 1] and validation_list[-4:] == ['tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'lF']:
                next_2_mere = error_segment[item[2] * 2:item[2] * 2 + 2]  # last,     current,    next,  next_next
                next_next_2_mere = error_segment[item[3] * 2:item[3] * 2 + 2]

                for variation in map_library[last_2_mere]:
                    if variation[0] == current_2_mere[0] or variation[1] == current_2_mere[1]:
                        for variation_next in map_library[variation]:
                            if variation_next[0] == next_2_mere[0] or variation_next[1] == next_2_mere[1]:
                                for variation_next_next in map_library[variation_next]:
                                    if variation_next_next[0] == next_next_2_mere[0] or variation_next_next[1] == next_next_2_mere[1]:
                                        possibilities.append(last_2_mere + variation + variation_next + variation_next_next)

                if possibilities != []:
                    trouble_parts.append([item, "Last 2-mere: 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'lF'", possibilities])


            ###############################
            ### 2-mere in Mitte falsch ####
            ###############################
            elif len(item) == 2:
                # last = 'tmT_nmF' --> richtig, current = 'tmF_nmT' --> falsch, danach wieder richtig --> alle möglichkeiten suchen
                for variation in map_library[last_2_mere]:
                    possibilities.append(last_2_mere + variation)

                # last = 'tmT_nmF' --> falsch, current = 'tmF_nmT' --> richtig, danach wieder richtig --> alle möglichkeiten suchen
                for i, items in enumerate(list_of_value):
                    if current_2_mere in items:
                        possibilities.append(list_of_key[i] + current_2_mere)

                # item anpassen nötig --> 1 davor falsch
                if possibilities != []:
                    trouble_parts.append([item, "2-mere wrong: 'tmT_nmF', 'tmF_nmT' ", possibilities])


            # T F T, 'tmT-nmF', 'tmF-nmF', 'tmF-nmT'
            elif len(item) == 3:
                next_2_mere = error_segment[item[2] * 2:item[2] * 2 + 2]

                for variation in map_library[last_2_mere]:
                    for j, rule_2_mere in enumerate(map_library[variation]):
                        # Rule in Variation stimmt überein → rule_2_mere als current_2_mere annehmen
                        if rule_2_mere == next_2_mere:  # variation ist der eigentliche 2-mere
                            # Variation aussortieren --> nur eine Substitution angenommen --> eine Base sollte passen!
                            if variation[0] == current_2_mere[0] or variation[1] == current_2_mere[1]:  # Level 1 diese Zeile wegnehmen
                                possibilities.append(last_2_mere + variation + next_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "2-mere wrong: 'tmT-nmF', 'tmF-nmF', 'tmF-nmT' ", possibilities])


            ################################
            ### 2 2-mere in Mitte falsch ###
            ################################
            # 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'     T F F T
            elif len(item) == 4:
                errorneous_part_1 = error_segment[item[1] * 2:item[1] * 2 + 2]
                errorneous_part_2 = error_segment[item[2] * 2:item[2] * 2 + 2]
                third_2_mere = error_segment[item[3] * 2:item[3] * 2 + 2]

                # Variationen 1-mere
                variation_option_1 = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in map_library[last_2_mere]:
                    if variation[0] == errorneous_part_1[0] or variation[1] == errorneous_part_1[1]:  # Level 1 diese Zeile wegnehmen
                        variation_option_1.append(variation)

                # Variationen 2-mere
                variation_option_2 = []
                for i, items in enumerate(list_of_value):
                    if third_2_mere in items:
                        if list_of_key[i][0] == errorneous_part_2[0] or list_of_key[i][1] == errorneous_part_2[1]:  # Level 1 diese Zeile wegnehmen
                            variation_option_2.append(list_of_key[i])

                # Variation 1 in Variation 2?
                possibilities = []
                for items in variation_option_1:
                    for j, rule_2_mere in enumerate(map_library[items]):
                        if rule_2_mere in variation_option_2:
                            possibilities.append(last_2_mere + items + rule_2_mere + third_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "2 2-mere wrong: 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'", possibilities])


            # 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'     T F T F T
            elif len(item) == 5:
                first_2_mere = error_segment[item[0] * 2:item[0] * 2 + 2]
                errorneous_part_1 = error_segment[item[1] * 2:item[1] * 2 + 2]
                middle_2_mere = error_segment[item[2] * 2:item[2] * 2 + 2]
                errorneous_part_2 = error_segment[item[3] * 2:item[3] * 2 + 2]
                last_2_mere = error_segment[item[4] * 2:item[4] * 2 + 2]

                # Möglichkeiten 1. fehlerhaftes 2-mere + passt es zum mittleren
                variation_option_1 = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in map_library[first_2_mere]:
                    if variation[0] == errorneous_part_1[0] and middle_2_mere in map_library[variation] or \
                            variation[1] == errorneous_part_1[1] and middle_2_mere in map_library[variation]:
                        variation_option_1.append(variation)

                # Möglichkeiten 2. fehlerhaftes 2-mere + passt es zum mittleren
                variation_option_2 = []
                for i, items in enumerate(list_of_value):
                    if last_2_mere in items:  # welcher key für last 2-mere
                        if list_of_key[i][0] == errorneous_part_2[0] or list_of_key[i][1] == errorneous_part_2[1]:  # eine Base sollte bei error-mere-2 passen
                            for j, mid in enumerate(list_of_value):
                                if list_of_key[i] in mid:  # welcher key für gefundene variation von error-mere-2
                                    if list_of_key[j] == middle_2_mere:  # muss mit mitte übereinstimmen
                                        variation_option_2.append(list_of_key[i])

                # Variation 1 in Variation 2?
                possibilities = []
                if variation_option_1 != [] and variation_option_2 != []:
                    for option_1 in variation_option_1:
                        for option_2 in variation_option_2:
                            possibilities.append(first_2_mere + option_1 + middle_2_mere + option_2 + last_2_mere)

                    trouble_parts.append([item, "2 2-mere wrong: 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'", possibilities])

        ###############################
        ## Sequenzen zusammensetzten ##
        ###############################
        sequence_parts = []
        start = 0
        if len(trouble_parts) == 0:
            decoded_sequences = []
        else:
            for trouble in trouble_parts:
                sequence_parts.append(error_segment[start:trouble[0][0] * 2])
                start = trouble[0][-1] * 2 + 2
            if start != len(error_segment):  # Einfügen des letzten Sequenzabschnittes, wenn nötig
                sequence_parts.append(error_segment[start:])

            possible_sequences = [sequence_parts[0]]
            # Fehler mit gefundenen Möglichkeiten ersetzen
            for i, trouble in enumerate(trouble_parts):  # trouble parts durchgehen + möglichkeiten zu possible_sequences hinzufügen und alte sequenz löschen
                new = []
                for seq in possible_sequences:
                    for items in trouble[2]:
                        if i + 1 < len(sequence_parts):
                            new.append(seq + items + sequence_parts[i + 1])
                        else:
                            new.append(seq + items)
                possible_sequences = new

            # Sequenzen übersetzen
            decoded_sequences = func.determine_correct_segments(possible_sequences, rs_codec_value, segment_count)

    else:
        decoded_sequences = []

    return decoded_sequences


def soft_level_1(error_segment, segment_count, rs_codec_value):
    """
    It is assumed that only substitutions are present. If inconsistencies are present, the mapping scheme is used to correct these with the available options.
    The correction is based on that of level 0, but without the restriction on the matching base required per 2-mere. In addition, the correction of 5
    consecutive incorrect 2-mers has been changed so that the 3 middle ones are all assumed to be incorrect.

    Args:
        error_segment: Input the errornous segment to correct.
        segment_count: Input the segment number of the segment for the DMR translation scheme.
        rs_codec_value: Input the used codec value for determination of the correct segments.

    Returns:
        List with the decoded sequence with were checked using the check_dmr_scheme.determine_correct_segments(possible_sequences, rs_codec_value,
        segment_count) function.
    """

    validation_list = func.validate_list(error_segment, segment_count)
    inconsistency_indices = func.find_inconsistency_indices_new(validation_list)
    inconsistency_count = len(inconsistency_indices)
    neighbouring_inconsistencies_list = func.get_inconsistency_neighbours_new(validation_list)

    for i, items in enumerate(neighbouring_inconsistencies_list):

        split_point = ['tmF_nmT']  # Definiere die Stellen, an denen die Liste geteilt werden soll
        result_list = []  # Initialisiere eine leere Liste für die geteilten Teillisten
        current_sublist = []  # Initialisiere eine Variable zum Zwischenspeichern der aktuellen Teilliste
        val = validation_list[items[0]:items[-1] + 1]
        for num, element in enumerate(val):
            # Füge das aktuelle Element zur aktuellen Teilliste hinzu
            current_sublist.append(items[num])

            # Überprüfe, ob das aktuelle Element eine der Split-Stellen ist
            if element in split_point:
                result_list.append(current_sublist)
                current_sublist = []

        if result_list != []:

            if current_sublist != []:
                result_list.append(current_sublist)

            neighbouring_inconsistencies_list.pop(i)
            for a in range(len(result_list), 0, -1):
                neighbouring_inconsistencies_list.insert(i, result_list[a - 1])

    len_neighbouring = [len(items) for items in neighbouring_inconsistencies_list]
    max_neighbour = max(len_neighbouring) if len_neighbouring else 0

    # keine Insertion/Deletion vorhanden, es müssen Fehler vom Mapping Schema da sein, höchstens 5 Inkonsistenzen hintereinander
    if inconsistency_count >= 1 and len(error_segment) % 2 == 0 and max_neighbour <= 5:

        #####################
        ### Vorbereitung ####
        #####################
        root_dir = os.path.abspath(os.curdir)
        with open(f"{root_dir}\\mapping_table_dmr.json") as f:
            data = json.load(f)
        map_library_start_2_mere = data["initial_2mer"]
        list_of_start_value = list(map_library_start_2_mere.values())
        map_library = data["map_library"]
        list_of_key = list(map_library.keys())
        list_of_value = list(map_library.values())
        trouble_parts = []

        # X
        for item in neighbouring_inconsistencies_list:
            last_2_mere = error_segment[item[0] * 2:item[0] * 2 + 2]
            current_2_mere = error_segment[item[1] * 2:item[1] * 2 + 2]
            possibilities = []

            ############################
            ### Start 2-mere falsch ####
            ############################
            # Start: F T
            if item == [0, 1] and validation_list[0: 2] == ['sF', 'tmF_nmT']:

                start_possibilities = list_of_start_value[segment_count % 4]
                for i, mere in enumerate(list_of_value):
                    if error_segment[2:4] in mere and list_of_key[i] in start_possibilities or \
                            error_segment[2:4] in mere and list_of_key[i] in start_possibilities:
                        possibilities.append(list_of_key[i] + error_segment[2:4])

                if possibilities != []:
                    trouble_parts.append([item, "First 2-mere: F, T", possibilities])

            # X
            # Start: F F T
            elif item == [0, 1, 2] and validation_list[0: 3] == ['sF', 'tmF_nmF', 'tmF_nmT']:
                next_2_mere = error_segment[item[2] * 2:item[2] * 2 + 2]

                # Variationen 1-mere
                variation_option_1 = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in list_of_start_value[segment_count % 4]:
                    variation_option_1.append(variation)

                # Variationen 2-mere
                variation_option_2 = []
                for i, items in enumerate(list_of_value):
                    if next_2_mere in items:
                        variation_option_2.append(list_of_key[i])

                # Variation 1 in Variation 2?
                possibilities = []
                for items in variation_option_1:
                    for j, rule_2_mere in enumerate(map_library[items]):
                        if rule_2_mere in variation_option_2:
                            possibilities.append(items + rule_2_mere + next_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "First 2-mere: F, F, T", possibilities])

            # X
            # 'sF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'
            elif item == [0, 1, 2, 3] and validation_list[0: 4] == ['sF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT']:

                errorneous_part_1 = error_segment[item[1] * 2:item[1] * 2 + 2]
                errorneous_part_2 = error_segment[item[2] * 2:item[2] * 2 + 2]
                third_2_mere = error_segment[item[3] * 2:item[3] * 2 + 2]

                # möglichen Start suchen
                start = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in list_of_start_value[segment_count % 4]:
                    start.append(variation)

                # Variationen 1-mere
                variation_option_1 = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for items in start:
                    for variation in map_library[items]:
                        variation_option_1.append(items + variation)

                # Variationen 2-mere
                variation_option_2 = []
                for i, items in enumerate(list_of_value):
                    if third_2_mere in items:
                        variation_option_2.append(list_of_key[i])

                # Variation 1 in Variation 2?
                possibilities = []
                for items in variation_option_1:
                    for j, rule_2_mere in enumerate(map_library[items[2:]]):
                        if rule_2_mere in variation_option_2:
                            possibilities.append(items + rule_2_mere + third_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "First 2-mere: F, F, F, T", possibilities])

            # X
            ##############################
            ### letztes 2-mere falsch ####
            ##############################
            # 'tmT_nmF', 'tF'
            elif item == [(len(error_segment) / 2) - 2, (len(error_segment) / 2) - 1] and validation_list[-2: -1] == ['tmT_nmF', 'lF']:
                for variation in map_library[last_2_mere]:
                    possibilities.append(last_2_mere + variation)

                if possibilities != []:
                    trouble_parts.append([item, "Last 2-mere: 'tmT_nmF', 'tF'", possibilities])

            # X
            # 'tmT_nmF', 'tmF_nmF', 'tF'
            elif item == [(len(error_segment) / 2) - 3, (len(error_segment) / 2) - 2, (len(error_segment) / 2) - 1] and \
                    validation_list[-3: -1] == ['tmT_nmF', 'tmF_nmF', 'tF']:
                next_2_mere = error_segment[item[2] * 2:item[2] * 2 + 2]

                # T F T
                for variation in map_library[last_2_mere]:
                    for j, rule_2_mere in enumerate(map_library[variation]):
                        if rule_2_mere == next_2_mere:  # variation ist der eigentliche 2-mere
                            possibilities.append(last_2_mere + variation + next_2_mere)

                # T F F   'tmT_nmF', 'tmF_nmF', 'lF'
                for variation in map_library[last_2_mere]:
                    for rule_2_mere in map_library[variation]:
                        possibilities.append(last_2_mere + variation + rule_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "Last 2-mere: 'tmT_nmF', 'tmF_nmF', 'tF'", possibilities])

            ###############################
            ### 2-mere in Mitte falsch ####
            ###############################

            # F, T  'tmT_nmF', 'tmF_nmT' --> 1. falsch passt jedoch ins Schema
            elif len(item) == 2:

                # last = 'tmT_nmF' --> richtig, current = 'tmF_nmT' --> falsch, danach wieder richtig --> alle möglichkeiten suchen
                for variation in map_library[last_2_mere]:
                    possibilities.append(last_2_mere + variation)

                # last = 'tmT_nmF' --> falsch, current = 'tmF_nmT' --> richtig, danach wieder richtig --> alle möglichkeiten suchen
                for i, items in enumerate(list_of_value):
                    if current_2_mere in items:
                        possibilities.append(list_of_key[i] + current_2_mere)

                # item anpassen nötig --> 1 davor falsch
                if possibilities != []:
                    trouble_parts.append([item, "2-mere wrong: 'tmT_nmF', 'tmF_nmT' ", possibilities])


            # T F T, 'tmT-nmF', 'tmF-nmF', 'tmF-nmT'
            elif len(item) == 3:
                next_2_mere = error_segment[item[2] * 2:item[2] * 2 + 2]

                for variation in map_library[last_2_mere]:
                    for j, rule_2_mere in enumerate(map_library[variation]):
                        # Rule in Variation stimmt überein → rule_2_mere als current_2_mere annehmen
                        if rule_2_mere == next_2_mere:  # variation ist der eigentliche 2-mere
                            possibilities.append(last_2_mere + variation + next_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "2-mere wrong: 'tmT-nmF', 'tmF-nmF', 'tmF-nmT' ", possibilities])


            ################################
            ### 2 2-mere in Mitte falsch ###
            ################################
            # 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'     T F F T
            elif len(item) == 4:
                errorneous_part_1 = error_segment[item[1] * 2:item[1] * 2 + 2]
                errorneous_part_2 = error_segment[item[2] * 2:item[2] * 2 + 2]
                third_2_mere = error_segment[item[3] * 2:item[3] * 2 + 2]

                # Variationen 1-mere
                variation_option_1 = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in map_library[last_2_mere]:
                    variation_option_1.append(variation)

                # Variationen 2-mere
                variation_option_2 = []
                for i, items in enumerate(list_of_value):
                    if third_2_mere in items:
                        variation_option_2.append(list_of_key[i])

                # Variation 1 in Variation 2?
                possibilities = []
                for items in variation_option_1:
                    for j, rule_2_mere in enumerate(map_library[items]):
                        if rule_2_mere in variation_option_2:
                            possibilities.append(last_2_mere + items + rule_2_mere + third_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "2 2-mere wrong", possibilities])

            # X
            # 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'     T F F F T
            elif len(item) == 5:
                first_2_mere = error_segment[item[0] * 2:item[0] * 2 + 2]
                last_2_mere = error_segment[item[4] * 2:item[4] * 2 + 2]

                # Möglichkeiten 1. fehlerhaftes 2-mere
                variation_option_1 = []  # Variation durchprobieren und schauen ob darauffolgender 2-mere im Schema ist
                for variation in map_library[first_2_mere]:
                    variation_option_1.append(variation)

                # Möglichkeiten 3. fehlerhaftes 2-mere
                variation_option_2 = []
                for i, items in enumerate(list_of_value):
                    if last_2_mere in items:  # welcher key für last 2-mere
                        variation_option_2.append(list_of_key[i])

                # Möglichkeiten 2. fehlerhaftes 2-mere
                possibilities = []
                for items in variation_option_1:
                    for i in map_library[items]:
                        for j in map_library[i]:
                            if j in variation_option_2:
                                possibilities.append(first_2_mere + items + i + j + last_2_mere)

                if possibilities != []:
                    trouble_parts.append([item, "2 2-mere wrong: 'tmT_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmF', 'tmF_nmT'", possibilities])

        ###############################
        ## Sequenzen zusammensetzten ##
        ###############################
        sequence_parts = []
        start = 0
        if len(trouble_parts) == 0:
            decoded_sequences = []
        else:
            for trouble in trouble_parts:
                sequence_parts.append(error_segment[start:trouble[0][0] * 2])
                start = trouble[0][-1] * 2 + 2
            if start != len(error_segment):  # Einfügen des letzten Sequenzabschnittes, wenn nötig
                sequence_parts.append(error_segment[start:])

            possible_sequences = [sequence_parts[0]]
            # Fehler mit gefundenen Möglichkeiten ersetzen
            for i, trouble in enumerate(
                    trouble_parts):  # trouble parts durchgehen und die möglichkeiten zu possible_sequences hinzufügen und alte sequenz löschen
                new = []
                for seq in possible_sequences:
                    for items in trouble[2]:
                        if i + 1 < len(sequence_parts):
                            new.append(seq + items + sequence_parts[i + 1])
                        else:
                            new.append(seq + items)
                possible_sequences = new

            # Sequenzen übersetzen
            decoded_sequences = func.determine_correct_segments(possible_sequences, rs_codec_value, segment_count)

    else:
        decoded_sequences = []

    return decoded_sequences


def soft_level_2(error_segment, segment_count, rs_codec_value):
    """
    It is assumed that only substitutions are present and that no correction with level 0 and 1 has worked. Here, the locations of the inconsistencies are
    first determined without separating them as in levels 0 and 1. Due to substitutions, errors could occur that fit into the mapping scheme and thus output
    false inconsistencies. In this case, a correction should be made by attempting to correct the list with the longest inconsistency. This involves going
    through the 2-mers and making individual substitutions by substituting each possible 2-mere (16 variants in total) at the position currently under
    consideration. The possibilities are corrected in level 0. If decodeable segments can be found, these are returned. If not, it continues with the next
    position. If no correctable segments are found with the substitution and level 0, the whole process is repeated with substitutions and a correction in
    level 1.

    Args:
        error_segment: Input the errornous segment to correct.
        segment_count: Input the segment number of the segment for the DMR translation scheme.
        rs_codec_value: Input the used codec value for determination of the correct segments.

    Returns:
        List with the decoded sequence with were checked using the check_dmr_scheme.determine_correct_segments(possible_sequences, rs_codec_value,
        segment_count) function.
    """

    validation_list = func.validate_list(error_segment, segment_count)
    inconsistency_indices = func.find_inconsistency_indices_new(validation_list)
    inconsistency_count = len(inconsistency_indices)
    neighbouring_inconsistencies_list = func.get_inconsistency_neighbours_new(validation_list)
    decoded_sequences = []

    if inconsistency_count >= 1 and len(error_segment) % 2 == 0:  # ansonsten Insertion/Deletion vorhanden, es müssen Fehler vom Mapping Schema da sein

        len_neighbouring = [len(items) for items in neighbouring_inconsistencies_list]
        if neighbouring_inconsistencies_list != []:
            max_neighbour = max(len_neighbouring) if len_neighbouring else 0
            index_neighbour = len_neighbouring.index(max_neighbour)

        # zuerst mit Level 0 probieren --> genauer
        for place in range(1, max_neighbour):
            option_list = \
                func.single_substitution_correction(error_segment, neighbouring_inconsistencies_list[index_neighbour][0] + place)[0]

            for option in option_list:
                decoded_sequences = soft_level_0(option, segment_count, rs_codec_value)
                if decoded_sequences != []:
                    return decoded_sequences

        # danach mit Level 1 probieren --> offener
        for place in range(1, max_neighbour):
            option_list = \
                func.single_substitution_correction(error_segment, neighbouring_inconsistencies_list[index_neighbour][0] + place)[0]

            for option in option_list:
                decoded_sequences = soft_level_1(option, segment_count, rs_codec_value)
                if decoded_sequences != []:
                    return decoded_sequences

    return decoded_sequences


def soft_level_3(error_segment, segment_count, rs_codec_value):
    """
    TO DO

    Args:
        error_segment: Input the errornous segment to correct.
        segment_count: Input the segment number of the segment for the DMR translation scheme.
        rs_codec_value: Input the used codec value for determination of the correct segments.

    Returns:
        List with the decoded sequence with were checked using the check_dmr_scheme.determine_correct_segments(possible_sequences, rs_codec_value,
        segment_count) function.
    """

    decoded_sequences = []
    return decoded_sequences





