import os
import sys
import optparse
import Levenshtein
import main_functions as do
import functions as func

parser = optparse.OptionParser()

parser.add_option("-o", "--out", metavar=" ", dest="out_path", default="output\\", help="Output filepath from project directory")
parser.add_option("-l", "--logo", metavar=" ", dest="use_logo", default="short", help="long or short", type=str)
parser.add_option("-p", "--percent", metavar=" ", dest="error_rates", default=0.01, help="List with error rates with 2 decimals. Input list in format -p X,Y,Z")
parser.add_option("-t", "--trial",  metavar=" ", dest="trials", default=10, help="Number of trials", type=int)
parser.add_option("-e", "--encoding", metavar=" ", dest="correction", default="DMR", help="Correction Mode: No, No_with_spacer, RS or RS_spacer, DMR, RS_segmented_packbits or DMR_segmented_packbits", type=str)
parser.add_option("-f", "--error", metavar=" ", dest="error_type", default="subs", help="all, subs, ins, del", type=str)
parser.add_option("-r", "--rscodec",  metavar=" ", dest="rs_codec_value", default=32, help="rs_codec_value", type=int)
parser.add_option("-c", "--mincodec",  metavar=" ", dest="min_codec_value", default=4, help="min_codec_value", type=int)
parser.add_option("-s", "--minseg",  metavar=" ", dest="min_segment_length", default=25, help="min_segment_length. Input list in format -s X,Y,Z")

if len(sys.argv) > 1:  # Start in the terminal when arguments are entered
    (options, args) = parser.parse_args()
    options.min_segment_length = [int(segment) for segment in options.min_segment_length.strip('[]').split(',')]    # Convert input to list
    print("min_segment_length", options.min_segment_length)
    options.error_rates = [float(segment) for segment in options.error_rates.strip('[]').split(',')]


else:   # If no arguments were passed from the command line, set the arguments here directly in the code
    options = optparse.Values()
    options.correction = "RS"
    options.out_path = "output\\"
    options.use_logo = "short"
    options.error_rates = [0.02]
    options.trials = 10
    options.error_type = "subs"
    options.rs_codec_value = 32
    options.min_codec_value = 1
    options.min_segment_length = [20, 25]    # [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]


#####################
##### Parameter #####
#####################
seed = 1

if options.correction == "RS" or options.correction == "RS_spacer":
    segment_length = 600
elif options.correction == "RS_segmented_packbits":
    segment_length = 800
elif options.correction == "DMR":
    segment_length = 64
elif options.correction == "DMR_segmented_packbits":
    segment_length = 120

rs_symbols_percent = 0.1254  # 12,54%  --> RS(255/223/32) --> 255 = 100% , 32 = 12,54%

if options.use_logo == "long":
    picture_path = "\\Logo_long__binary_.png"
    picture_width = 1280
    picture_height = 210
    bitlen = 16 if options.correction == "RS" else 8
elif options.use_logo == "short":
    picture_path = "\\Logo_short__binary_.png"
    picture_width = 47
    picture_height = 47
    bitlen = 10 if options.correction == "RS" else 8


#####################
##### Encoding ######
#####################
for length in options.min_segment_length:

    if options.correction == "RS_segmented_packbits" or options.correction == "DMR_segmented_packbits":
        dna_message_list, binary_segment_list = \
            do.encoding_picture_segmented_packbits(seed, bitlen, segment_length, rs_symbols_percent, options.correction,
                                                   logging_filepath=options.out_path + options.correction + "//",
                                                   logging_filename="length_" + str(length) + "_original_dna_encoding.log", logging_num_elements=100,
                                                   logging_num_segments=10,
                                                   rs_codec_value=options.rs_codec_value, min_codec_value=options.min_codec_value, min_segment_length=length,
                                                   picture_path=picture_path, picture_type="binary", show_image=False, verbose=True)
    else:
        dna_message_list, binary_segment_list = do.encoding_picture(seed, bitlen, segment_length, rs_symbols_percent, options.correction, False,
                                                                    logging_filepath=options.out_path + options.correction + "//",
                                                                    logging_filename="length_" + str(length) + "_original_dna_encoding.log",
                                                                    logging_num_elements=100, logging_num_segments=10,
                                                                    rs_codec_value=options.rs_codec_value, min_codec_value=options.min_codec_value,
                                                                    min_segment_length=length, picture_path=picture_path, picture_type="binary",
                                                                    random_mask=True, show_image=False, verbose=True, do_xor=False)

    # Calculate length and output with hint + spacer
    dna = "XXXXXX".join(dna_message_list)
    dna_length = len(dna) - 6 * (len(dna_message_list) - 1)

    with open(options.out_path + "//" + options.correction + "//" + "length_" + str(length) + "_original_dna.fasta", "w") as outfile:
        if options.correction == "RS" or options.correction == "No":
            outfile.write("Länge der erzeugten DNA: " + str(dna_length) + "\n")
        else:
            outfile.write("Länge der erzeugten DNA: " + str(dna_length) + " + " + str(len(dna_message_list)-1) + " Spacersequenzen.\n")
        outfile.write("Ohne Spacer einzubeziehen ergibt sich eine Nettoinformationsdichte von = " + str((picture_width * picture_height)/dna_length) + "\n")
        outfile.write("\n".join([dna[n:n + 60] for n in range(0, len(dna), 60)]))


    ###################################
    ##### Add Errors and Decoding #####
    ###################################
    for error in options.error_rates:
        for trial in range(options.trials):

            loop = True
            while loop:

                try:
                    ############################
                    ##### Error simulation #####
                    ############################
                    print("Adding errors to sequence")
                    dna_mut = ""
                    if error == 0.00:
                        dna_mut = dna
                    else:
                        if options.error_type == "all":
                            # dna_mut = err.binom_mutations(dna[0], subs, ins, dels, binom=False)
                            dna_mut = func.binom_mutations_with_spacer_ignorance(dna, dna_length, error / 3, error / 3, error / 3, binom=False)
                        elif options.error_type == "subs":
                            dna_mut = func.binom_mutations_with_spacer_ignorance(dna, dna_length, error, 0, 0, binom=False)
                        elif options.error_type == "ins":
                            dna_mut = func.binom_mutations_with_spacer_ignorance(dna, dna_length, 0, error, 0, binom=False)
                        elif options.error_type == "del":
                            dna_mut = func.binom_mutations_with_spacer_ignorance(dna, dna_length, 0, 0, error, binom=False)


                    # Split dna_mut by X so that sections can be decoded
                    dna_mut_list = dna_mut.split("X")
                    dna_mut_list = [item for item in dna_mut_list if item != ""]

                    # Calculate length and output with hint + spacer
                    dna_mut_length = sum(len(item) for item in dna_mut_list)

                    with open(options.out_path + options.correction + "//" + "length_" + str(length) + "_" + options.error_type +
                              "_error-rate_" + str(error) + "_mut_dna.fasta", "a") as outfile:
                        if options.correction == "RS" or options.correction == "No":
                            outfile.write("Länge der mutierten DNA: " + str(dna_mut_length) + "\n")
                        else:
                            outfile.write("Länge der mutierten DNA: " + str(dna_mut_length) + " + " + str(len(dna_message_list) - 1) + " Spacersequenzen.\n")

                        outfile.write("%s\n" % dna_mut_list)


                    ####################
                    ##### Decoding #####
                    ####################
                    if options.correction == "RS_segmented_packbits" or options.correction == "DMR_segmented_packbits":
                        bitseq, picture, _, _ = \
                            do.picture_decode_segmented_packbits(dna_mut_list, seed, bitlen, rs_symbols_percent, segment_length, options.correction, True,
                                                                 picture_width, picture_height, options.out_path + options.correction + "//",
                                                                 "length_" + str(length) + "_" + options.error_type +
                                                                 "_error-rate_" + str(error) + "_trial_" + str(trial) + "_decoded_picture",
                                                                 options.out_path + options.correction + "//",
                                                                 "length_" + str(length) + "_" + options.error_type + "_error-rate_" + str(error) + ".log",
                                                                 logging_num_elements=100, logging_num_segments=10, picture_mode="binary",
                                                                 correction=True, verbose=True)

                    else:
                        bitseq, picture = do.picture_decode(dna_mut_list, seed, bitlen, rs_symbols_percent, segment_length, options.correction, False,
                                                            picture_width, picture_height, options.out_path + options.correction + "//",
                                                            "length_" + str(length) + "_" + options.error_type + "_error-rate_" + str(error) +
                                                            "_trial_" + str(trial) + "_decoded_picture", options.out_path + options.correction + "//",
                                                            "length_" + str(length) + "_" + options.error_type + "_error-rate_" + str(error) + ".log",
                                                            logging_num_elements=100, logging_num_segments=10, rs_codec_value=options.rs_codec_value,
                                                            min_codec_value=options.min_codec_value, min_segment_length=length,
                                                            picture_mode="binary", random_mask=True, verbose=True, do_xor=False,
                                                            correction=True)

                    loop = False
                except:
                    loop = True


            with open(options.out_path + options.correction + "//" + "length_" + str(length) + "_" + options.error_type +
                      "_error-rate_" + str(error) + "_bitseq.fasta", "a") as outfile:
                outfile.write("Länge der dekodierten Binärsequenz: " + str(len(bitseq)) + "\n")
                outfile.write("\n".join([bitseq[n:n + picture_width] for n in range(0, len(bitseq), picture_width)]))
                outfile.write("\n\n")

            ########################
            ##### Edit distance ####
            ########################
            print("Start Edit Distance calculation: assembled Seq")
            # load ref pictures binary
            ROOT_DIR = os.path.abspath(os.curdir)
            pic_ref = func.read_1bit_image(ROOT_DIR + picture_path, show_image=False)
            pic_ref = "".join([str(bit) for bit in pic_ref])
            pic_dec = func.read_1bit_image(ROOT_DIR + "//" + options.out_path + options.correction + "//" + "length_" + str(length) + "_" +
                                           options.error_type + "_error-rate_" + str(error) + "_trial_" + str(trial) + "_decoded_picture.png", show_image=False)
            pic_dec = "".join([str(bit) for bit in pic_dec])
            # calculation
            edit_scores_dna = func.score_pairs_fast([(dna, dna_mut, "Encoded & mutated DNA:")], as_percent=True)
            edit_scores_bit = func.score_pairs_fast([(pic_ref, pic_dec, "Encoded & decoded picture")])[0] if pic_ref != bitseq else ("Encoded & decoded picture", 1.0)

            # Edit Add distance and name to file
            with open(options.out_path + options.correction + "//" + "length_" + str(length) + "_" + options.error_type
                      + "_error-rate_" + str(error) + "_edit_distances.txt", "a") as file:
                file.write("\nTrial: " + str(trial) + ", " + str(edit_scores_dna[0]) + ", " + str(edit_scores_bit))

