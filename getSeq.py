#!/usr/bin/env python3
#  -*- coding: utf-8 -*-

# 2017-03-23
# STB

import os
import sys
import pandas as pd
import regex as re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from file_read_backwards import FileReadBackwards

prog_version = "0.11"

description = """Extract DNA sequences from larger contexts like genome
sequences in the form of a FASTA file. The user can supply a
GFF3 file containing the FASTA file's annotation in order to then extract
DNA sequences by their gene symbols (as found in the GFF file). It is
possible to specify which annotated element of a given gene shall be
extracted (e.g. gene, mRNA, CDS).
Alternatively, the user can manually set start and end positions of a
sequence and its strandedness in order to retrieve that sequence from the
FASTA file.
In both usage scenarios it is possible to define a number of nucleotides to be
extracted up- and/or downstream of the specified annotated element. 
Additionally, it is also possible to just extract up- or downstream 
sequences without the related annotated element's sequence itself 
(e.g. in order to create promoteromes)."""

MIT_license = """Copyright 2017 Sven T. Bitters (sven.bitters@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


def parse_input():
    parser = ArgumentParser(description=description,
                            usage="getSeq.py [-h] [-l]\n"
                                  "                 -g GENOME_FASTA -a GENOME_GFF\n"
                                  "                 [-t TARGET | -L TARGET_LIST] [-y TARGET_TYPE]\n"
                                  "                 [-m TARGET_TRANSCRIPT] [-u UPSTREAM_NUM] [-d DOWNSTREAM_NUM]\n"
                                  "                 [-f FORMAT] [--onlyupstream [-i INSIDE] | --onlydownstream]\n"
                                  "                 [-x EXCLUDE_REGEX] [-o OUTPUT_FILE]\n"
                                  "                 or\n"
                                  "                 -n FASTA_ID -c CHROMOSOME -s START -e END -r STRAND\n"
                                  "                 [-u UPSTREAM_NUM] [-d DOWNSTREAM_NUM] [-f FORMAT]\n"
                                  "                 [--onlyupstream [-i INSIDE] | --onlydownstream]\n"
                                  "                 [-o OUTPUT_FILE]")
    parser.add_argument("-l", "--license", action="store_true", default=False, dest="license",
                        help="show license and exit")
    parser.add_argument("-g", "--genome", type=str, dest="genome_fasta",
                        help="String. Mandatory. Path to a (multi-)FASTA file containing DNA sequences.")
    parser.add_argument("-a", "--annotation", type=str, dest="genome_gff",
                        help="String. Mandatory. Path to a GFF3 file containing annotations corresponding to the FASTA file. "
                             "The GFF file must follow GFF3 file format specifications (http://www.sequenceontology.org/gff3.shtml) "
                             "with IDs following this pattern: 'ID=gene1337' or 'ID=cds80085' or 'ID=rna42' or 'ID=id31415'.")
    parser.add_argument("-o", "--output", type=str, dest="output_file", default="~/Documents/getSeq/",
                        help="String. Either the path to a directory where getSeq's output file will be saved or the complete path "
                             "for a file ending in '.fa'. If you do only supply a directory for output, getSeq will choose a file name "
                             "based on your query (e.g. if you supply getSeq with a list of targets via -L, the list's name will be "
                             "used in the file name). (default: ~/Documents/getSeq/<context>.fa)")
    group_input = parser.add_mutually_exclusive_group()
    group_input.add_argument("-t", "--target", default=False, dest="target",
                        help="String. ID of one gene as written in the the 'gene=' element in the attribute column of the annotation file. "
                             "This gene's sequence will be extracted.")
    group_input.add_argument("-L", "--list", default=False, dest="target_list",
                        help="String. Path to a file containing either gene IDs as written in the annotation file or specified regions in the "
                             "genome that will be extracted one after the other. The output DNA sequences will be saved to a multi-FASTA file.")
    parser.add_argument("-y", "--type", type=str, nargs="+", default=["gene", "single"], dest="target_type",
                        help="String. Annotation type as available in the genome's annotation file, e.g. 'gene', 'mRNA', 'exon'. "
                             "The DNA sequence to be extracted is based on the information in the specified gene ID's <type> annotation. "
                             "The way this option is executed is somewhat dependent on which other arguments you have supplied: "
                             "a) If you have specified a target gene via -t, only the specified <type> of the selected gene is extracted by getSeq; "
                             "b) If you have specified a list of target genes via -L, the <type> annotation of these genes will be extracted and "
                             "written in a multi-FASTA file; c) If you have neither specified a target gene via -t or a list of target genes via -L "
                             "ALL <type> of ALL genes in the annotation will be extracted. If in this case <type> does not exist for a given gene, "
                             "getSeq will revert to using the 'gene' annotation of that gene. In default mode ('single'), all <type> of each gene will be written "
                             "to a multi-FASTA file. Thus, if you have typed e.g. '-y CDS', all individual CDS annotation's DNA sequences will be output. "
                             "However, it is also possible to supply getSeq's -y option with a second argument to change this behavior. By typing "
                             "e.g. '-y mRNA span' you will receive one long DNA sequence spanning from the first to the last mRNA annotation for a "
                             "given gene - i.e. the most 5' 'start' annotation and the most 3' 'end' annotation of <type> will be used to generate this long "
                             "DNA sequence. Additionally, by typing e.g. '-y CDS merged' you can merge <type> annotations which share the same ID. "
                             "In effect, the DNA sequences corresponding to the individual <type> annotations will be concatenated in order to generate "
                             "one long continuous sequence. This option is only available for annotation elements sharing their ID - which is mostly the case with "
                             "'CDS' annotations because often CDS annotations that make up one transcript's ORF share IDs (i.e. all CDS annotations on "
                             "transcript 'rna42' share the ID 'cds1337'. (default: gene single)")
    parser.add_argument("-m", "--element", default=0, dest="target_transcript",
                        help="String or integer. If you have specified anything other than 'gene' via -y you can select which "
                             "transcript/CDS/exon you want to have returned - either by specifying its specific ID (e.g. rna42) or by its "
                             "number (e.g. 2 for the second exon). Annotations specified via -p must have individual IDs - if this is not the case "
                             "(like for CDS annotations) -m will select the set of all annotations sharing an ID. I.e. if you type '-y CDS.merged -p 2' "
                             "you will receive the merged DNA sequences of the second CDS annotation set for the given gene. If the specified element is not "
                             "present in your annotation, the next available element is used, i.e. if you have requested the 5th transcript of a "
                             "gene which gives rise to only 3 transcripts, you will receive the 3rd transcript. If you do not supply a target "
                             "gene via -t ALL annotations of the specified type will be considered by getSeq. The final output then depends on the mode "
                             "selected via -y X.<mode>. (default: 0 - this means 'all')")
    parser.add_argument("-u", "--upstream", default=0, type=int, dest="upstream_num",
                        help="Integer. Specify a number of nucleotides upstream of the annotated element to be included in the output. "
                             "If you are submitting a -y or -L request the specified number of nucleotides will be added to each annotated "
                             "element. (default: 0)")
    parser.add_argument("-d", "--downstream", default=0, type=int, dest="downstream_num",
                        help="Integer. Specify a number of nucleotides downstream of the annotated element to be included in the output. "
                             "If you are submitting a -y or -L request the specified number of nucleotides will be added to each annotated "
                             "element. (default: 0)")
    parser.add_argument("-f", "--format", default="fasta", type=str, dest="fa_format",
                        help="String. Specify in which format the output DNA sequences should be written to a file. You can either select to "
                             "write the DNA sequences in FASTA format (i.e. line-breaks every 60 characters) or in single-line format (i.e. "
                             "there will be a FASTA header in one line and the whole DNA sequence in the next). "
                             "Options: fasta, single-line. (default: fasta)")
    group_only = parser.add_mutually_exclusive_group()
    group_only.add_argument("--onlyupstream", action="store_true", default=False, dest="upstream_only",
                        help="Flag. If this flag is set, only the upstream DNA sequence specified via -u will be written to the output FASTA file "
                             "- thus, the annotated element's own DNA sequence will be ommited. Useful when creating promoteromes.")
    group_only.add_argument("--onlydownstream", action="store_true", default=False, dest="downstream_only",
                        help="Flag. If this flag is set, only the downstream DNA sequence specified via -d will be written to the output FASTA file "
                             "- thus, the annotated element's own DNA sequence will be ommited.")
    parser.add_argument("-i", "--inside", default=0, type=int, dest="inside",
                        help="Integer. If you have decided to return only the upstream sequences of your annotated elements via --onlyupstream "
                             "you can also specify a number of nucleotides to include of said annotated element in your output FASTA file. Thus, "
                             "if you for example type '-i 100' the first 100 nt of the annotated element(s) will be included in the FASTA file. (default: 0)")
    parser.add_argument("-x", "--exclude", default="", type=str, dest="exclude_regex",
                        help="String. Specify a Python-style regular expression in order to exclude genes with matching names from being processed.")
    parser.add_argument("-n", "--name", default="userdef_name.userdef_id", type=str, dest="fasta_id",
                        help="String. Manual sequence retrieval. Name to use in the output FASTA file's FASTA header.")
    parser.add_argument("-c", "--chromosome", default=False, dest="chromosome",
                        help="String. Manual sequence retrieval. Chromosome on which the sequence to be retrieved is located.")
    parser.add_argument("-s", "--start", default=-1, type=int, dest="start",
                        help="Integer. Manual sequence retrieval. Start position on the chromosome of the sequence to be retrieved.")
    parser.add_argument("-e", "--end", default=-1, type=int, dest="end",
                        help="Integer. Manual sequence retrieval. End position on the chromosome of the sequence to be retrieved.")
    parser.add_argument("-r", "--strand", default=False, type=str, dest="strand",
                        help="String. Manual sequence retrieval. Strand of the chromosome on which the sequence to be retrieved is located." )

    parserargs = parser.parse_args()

    try:
        if parserargs.license:
            print(MIT_license)

        else:
            input_fasta = parserargs.genome_fasta
            input_gff = parserargs.genome_gff

            if input_fasta is None and input_gff is None:
                parser.print_help()
                raise SystemExit

            else:
                output_path = parserargs.output_file
                target_list = parserargs.target_list
                target_element = parserargs.target
                target_type = parserargs.target_type
                target_transcript = parserargs.target_transcript
                with_upstream = parserargs.upstream_num
                with_downstream = parserargs.downstream_num
                fasta_id = parserargs.fasta_id
                mychromosome = parserargs.chromosome
                mystart = parserargs.start
                myend = parserargs.end
                mystrand = parserargs.strand
                dnaformat = parserargs.fa_format
                with_inside = parserargs.inside
                exclude_raw = parserargs.exclude_regex
                upstream_only = parserargs.upstream_only
                downstream_only = parserargs.downstream_only

                check_inputs(input_fasta, input_gff, output_path, target_list,
                             target_element, target_type, target_transcript,
                             with_upstream, with_downstream, upstream_only,
                             downstream_only, fasta_id, mychromosome, mystart,
                             myend, mystrand, dnaformat, with_inside,
                             exclude_raw)

                return input_fasta, input_gff, output_path, target_list, target_element, target_type, target_transcript, \
                       with_upstream, with_downstream, upstream_only, downstream_only, fasta_id, mychromosome, mystart, \
                       myend, mystrand, dnaformat, with_inside, exclude_raw

    except SystemExit:
        sys.exit()
    except TypeError:
        print("Error: Please supply a FASTA AND a GFF3 file.")
        sys.exit()


def check_inputs(input_fasta, input_gff, output_path, target_list, target_element, target_type, target_transcript,
                       with_upstream, with_downstream, upstream_only, downstream_only, fasta_id, mychromosome, mystart,
                       myend, mystrand, dnaformat, with_inside, exclude_raw):

    fa_test = 0
    with open(input_fasta, "r") as test_fa_handle:
        legal_fa_symbols = ["A", "T", "G", "C", "N", "U", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "-"]
        for line in test_fa_handle:
            if fa_test == 0:
                if line[0] != ">" or line[1] == " " or len(line.rstrip()) < 2:
                    print("Error: Input FASTA file does not conform to FASTA file type specifications: bad FASTA header")
                    raise SystemExit
            elif fa_test == 1:
                if line[0] != ";":
                    if not all(nuc in legal_fa_symbols for nuc in set(line.rstrip().upper())):
                        print("Error: Input FASTA file does not conform to FASTA file type specifications: bad nucleotides.")
                        raise SystemExit
                else:
                    fa_test -= 1
            else:
                break
            fa_test += 1

        if fa_test < 2:
            print("Error: Input FASTA file does not conform to FASTA file type specifications: no DNA sequence.")
            raise SystemExit

    gff_test = 0
    with open(input_gff, "r") as test_gff_handle:
        for line in test_gff_handle:
            if not line.startswith("#"):
                gff_line = line.split("\t")

                if len(gff_line) != 9 \
                   or gff_line[6] not in ["+", "-", ".", "?"] \
                   or "ID=" not in gff_line[8]:
                    print("Error: Input GFF3 file does not conform to GFF3 file type specifications: bad format.")
                    raise SystemExit

                gff_test += 1
                break

        if gff_test != 1:
            print("Error: Input GFF3 file does not conform to GFF3 file type specifications: file does not contain information.")
            raise SystemExit

    if input_gff is None and (target_element is not None or target_list is not None):
        print("Error: Please supply a GFF3 file.")
        raise SystemExit

    if (target_element is not None or target_list is not None) and \
       (mystart > -1 or myend > -1):
        print("Error: You cannot request a specific target to be extracted and at the same time provide a list of targets. "
              "Please remove one.")
        raise SystemExit

    if upstream_only and with_upstream == 0:
        print("Error: You have requested to only output upstream sequences but did not specify how long the upstream sequence should be. "
              "Please choose an upstream sequence length.")
        raise SystemExit
    elif downstream_only and with_downstream == 0:
        print("Error: You have requested to only output downstream sequences but did not specify how long the downstream sequence should be. "
              "Please choose a downstream sequence length.")
        raise SystemExit

    if (upstream_only is False and with_inside != 0) or \
       (downstream_only is True and with_inside != 0):
        print("Error: You cannot have an 'inside' sequence if you do not set the --upstreamonly flag.")
        raise SystemExit

    if with_upstream < 0:
        print("Error: Upstream sequence length cannot be < 0.")
        raise SystemExit
    elif with_downstream < 0:
        print("Error: Upstream sequence length cannot be < 0.")
        raise SystemExit
    elif with_inside < 0:
        print("Error: 'Inside' sequence length cannot be < 0.")
        raise SystemExit

    if dnaformat not in ["fasta", "single-line"]:
        print("Error: You have requested an unknown output format. Legal options: fasta, single-line")
        raise SystemExit

    if (mystart > -1 and myend == -1) or (mystart == -1 and myend > -1):
        print("Error: You are requesting a manually selected sequence but did not specify both, start and end value for your sequence. "
              "Please specify both, start and end value.")
        raise SystemExit

    if mystrand and mystrand not in ["-", "+"]:
        print("Error: You are requesting a manually selected sequence but specified an illegal value for strandedness. "
              "Please specify whether your sequence is on + or -.")
        raise SystemExit

    if mystart > -1 and myend > -1 and mychromosome and not mystrand:
        print("Error: You are requesting a manually selected sequence but did not specify your sequence's strandedness. "
              "Please specify whether your sequence is on + or -.")
        raise SystemExit
    elif mystart > -1 and myend > -1 and mystrand and not mychromosome:
        print("Error: You are requesting a manually selected sequence but did not specify on which chromosome your sequence is located. "
              "Please specify the chromosome's name.")
        raise SystemExit
    elif mystart == -1 and myend == -1 and mystrand and mychromosome:
        print("Error: You are requesting a manually selected sequence but did not specify start and end value for your sequence. "
              "Please specify start and end value.")
        raise SystemExit

    if len(target_type) == 1:
        target_type.append("single")

    if target_type[1] not in ["single", "span", "merged"]:
        print("Error: You have requested an unknown target mode. Legal options: single, span, merged")
        raise SystemExit


def read_and_mod_gff(annot_gff):
    # Read the GFF file and replace the values "attributes" column with just the IDs of the annotated elements
    # i.e. something like "ID=id364475;Parent=gene41724;Dbxref=GeneID:19989172;exon_number=1;gbkey=exon;gene=cox2"
    # becomes "id364475"

    id_regex = re.compile(r"(?<=ID=).+?(?=($|;))")
    parent_regex = re.compile(r"(?<=Parent=).+?(?=($|;))")
    loc_regex = re.compile(r"(?<=gene=).+?(?=($|;))")
    name_regex = re.compile(r"(?<=Name=).+?(?=($|;))")

    id_dict = {}
    last_name = ""

    gff_line_list = list()
    print("Reading GFF...")
    gff_in = FileReadBackwards(annot_gff, encoding="utf-8")
    for line in gff_in:

        if re.match("\w", line) and not line.startswith('#'):
            tab_elements = line.split("\t")
            type = tab_elements[2]
            attributes = tab_elements[-1]

            try:
                element_id = re.search(id_regex, attributes)
                element_id = element_id.group()
            except:
                element_id = "."

            if element_id == ".":
                id_dict[type] = id_dict.get(type, 0) + 1
                element_id = type + str(id_dict[type]-1)

            try:
                loc_name = re.search(loc_regex, attributes)
                loc_name = loc_name.group()
            except:
                loc_name = "."

            try:
                gene_name = re.search(name_regex, attributes)
                gene_name = gene_name.group()
                last_name = gene_name
            except:
                gene_name = "."

            if gene_name == ".":
                gene_name = last_name

            if gene_name != "." and loc_name == ".":
                loc_name = gene_name

            try:
                parent_gene = re.search(parent_regex, tab_elements[-1])
                parent_gene = parent_gene.group()
            except:
                parent_gene = "."

            if parent_gene == "." and gene_name != ".":
                parent_gene = gene_name

            tab_elements[-1] = element_id
            tab_elements += [parent_gene, loc_name]

            gff_line_list.append(tab_elements)

    gff_cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "id", "parent", "gene_name"]
    gff_df = pd.DataFrame(gff_line_list, columns=gff_cols)

    return gff_df


def find_target(input_fasta, gff_df, with_upstream, with_downstream, upstream_only, downstream_only, target_list,
                target_element, target_type, target_transcript, output_fasta, fasta_id, mychromosome, mystart, myend,
                mystrand, dnaformat, with_inside, exclude_regex):

    print("Finding requested annotation(s)...")

    gene_occurence_dict = dict()
    id_occurence_dict = dict()

    element_id = "id"

    target_mode = target_type[1]
    target_type = target_type[0]

    gff_df_gene = gff_df[(gff_df["type"] == "gene")]
    gff_df_subset = gff_df[(gff_df["type"] == target_type)]

    print("Writing FASTA...")

    with open(output_fasta, 'w') as output_fasta_handle:

        if target_element:
            chrom_dict = False

            element_id, gene_occurence_dict, id_occurence_dict = get_gff_target(input_fasta, gff_df_gene, gff_df_subset,
                                                                                target_element, target_transcript,
                                                                                chrom_dict, target_mode, with_upstream,
                                                                                with_downstream, output_fasta_handle,
                                                                                target_type, upstream_only,
                                                                                downstream_only, dnaformat, with_inside,
                                                                                exclude_regex, gene_occurence_dict,
                                                                                id_occurence_dict)


        elif mystart != -1:
            chrom_id = mychromosome
            startpos = int(mystart)
            endpos = int(myend)
            strandedness = str(mystrand)
            upstream = int(with_upstream)
            downstream = int(with_downstream)

            chrom_dict = False
            element_id = "userdef_id"
            element_name = "userdef_name"

            gene_occurence_dict_null, id_occurence_dict_null = extract_seq(input_fasta, chrom_dict, chrom_id, startpos,
                                                                           endpos, strandedness, element_name,
                                                                           element_id, upstream, downstream,
                                                                           gene_occurence_dict, id_occurence_dict,
                                                                           output_fasta_handle, target_type,
                                                                           upstream_only, downstream_only, dnaformat,
                                                                           with_inside, exclude_regex, fasta_id)

        elif target_list or target_type:
            chrom_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, format="fasta"))

            if target_list:
                with open(target_list, "r") as target_list_handle:
                    for line in target_list_handle:
                        line = line.rstrip()
                        line_list = line.split("\t")

                        if len(line_list) == 5 or len(line_list) == 7:
                            element_name = "."
                            element_id = str(line_list[0])
                            chrom_id = str(line_list[1])
                            startpos = int(line_list[2])
                            endpos = int(line_list[3])
                            strandedness = str(line_list[4])

                            if len(line_list) == 5:
                                upstream = int(with_upstream)
                                downstream = int(with_downstream)
                            elif len(line_list) == 7:
                                upstream = int(line_list[5])
                                downstream = int(line_list[6])

                            gene_occurence_dict, id_occurence_dict = extract_seq(input_fasta, chrom_dict, chrom_id,
                                                                                 startpos, endpos, strandedness,
                                                                                 element_name, element_id, upstream,
                                                                                 downstream, gene_occurence_dict,
                                                                                 id_occurence_dict,
                                                                                 output_fasta_handle, target_type,
                                                                                 upstream_only, downstream_only,
                                                                                 dnaformat, with_inside, exclude_regex)

                        else:
                            target_element = line_list[0]

                            if len(line_list) == 2 or len(line_list) == 4:
                                target_transcript = line_list[1]

                            if len(line_list) == 3:
                                upstream = int(line_list[1])
                                downstream = int(line_list[2])
                            elif len(line_list) == 4:
                                upstream = int(line_list[2])
                                downstream = int(line_list[3])
                            else:
                                upstream = int(with_upstream)
                                downstream = int(with_downstream)

                            element_id, gene_occurence_dict, id_occurence_dict = get_gff_target(input_fasta,
                                                                                                gff_df_gene,
                                                                                                gff_df_subset,
                                                                                                target_element,
                                                                                                target_transcript,
                                                                                                chrom_dict, target_mode,
                                                                                                upstream, downstream,
                                                                                                output_fasta_handle,
                                                                                                target_type,
                                                                                                upstream_only,
                                                                                                downstream_only,
                                                                                                dnaformat, with_inside,
                                                                                                exclude_regex,
                                                                                                gene_occurence_dict,
                                                                                                id_occurence_dict)


            elif target_type:
                for index, gff_line_gene in gff_df_gene.iterrows():

                    gene_target = gff_line_gene["gene_name"]

                    element_id, gene_occurence_dict, id_occurence_dict = get_gff_target(input_fasta, gff_df_gene,
                                                                                        gff_df_subset, gene_target,
                                                                                        target_transcript, chrom_dict,
                                                                                        target_mode, with_upstream,
                                                                                        with_downstream,
                                                                                        output_fasta_handle,
                                                                                        target_type, upstream_only,
                                                                                        downstream_only, dnaformat,
                                                                                        with_inside, exclude_regex,
                                                                                        gene_occurence_dict,
                                                                                        id_occurence_dict)

    return element_id


def get_gff_target(input_fasta, gene_gff_df, subset_gff_df, gene_input, target_transcript, chrom_dict, target_mode,
                   with_upstream, with_downstream, output_fasta_handle, target_type, upstream_only, downstream_only,
                   dnaformat, with_inside, exclude_regex, gene_occurence_dict, id_occurence_dict):

    gene_gff_line = gene_gff_df[gene_gff_df["gene_name"] == gene_input]
    subset_gff_lines = subset_gff_df[subset_gff_df["gene_name"] == gene_input]

    if subset_gff_lines.empty:
        subset_gff_lines = gene_gff_line
        target_type = "gene"
        target_mode = "single"

    subset_ids = sorted(set(subset_gff_lines["id"]))

    if target_transcript != -1 and target_type != "gene":
        try:
            if target_type == "mRNA":
                subset_gff_lines = subset_gff_lines.iloc[target_transcript]
            else:
                subset_ids = [subset_ids[target_transcript]]
        except IndexError:
            for ii in range(target_transcript, -1, -1):
                try:
                    if target_type == "mRNA":
                        subset_gff_lines = subset_gff_lines.iloc[ii]
                    else:
                        subset_ids = subset_ids[ii]
                    break
                except IndexError:
                    pass

        except TypeError:
            if target_type == "mRNA":
                subset_gff_lines = subset_gff_lines[subset_gff_lines["id"] == target_transcript]
            else:
                subset_ids = [subset_ids(subset_ids.index(target_transcript))]

    if not isinstance(subset_gff_lines, pd.DataFrame):
        subset_gff_lines = pd.DataFrame([subset_gff_lines])

    if target_mode == "single":
        for index, subset_line in subset_gff_lines.iterrows():
            element_name = str(subset_line["gene_name"])
            chrom_id = str(subset_line["seqid"])
            strandedness = str(subset_line["strand"])
            element_id = str(subset_line["id"])
            upstream = int(with_upstream)
            downstream = int(with_downstream)
            startpos = int(subset_line["start"])
            endpos = int(subset_line["end"])

            gene_occurence_dict, id_occurence_dict = extract_seq(input_fasta, chrom_dict, chrom_id, startpos, endpos,
                                                                 strandedness, element_name, element_id,
                                                                 upstream, downstream, gene_occurence_dict,
                                                                 id_occurence_dict, output_fasta_handle, target_type,
                                                                 upstream_only, downstream_only, dnaformat, with_inside,
                                                                 exclude_regex)

    elif target_mode == "span":
        start_list = list()
        end_list = list()

        if len(subset_gff_lines["start"]) > 1:
            for index, sub_line in subset_gff_lines.iterrows():
                start_list.append(int(sub_line["start"]))
                end_list.append(int(sub_line["end"]))
        else:
            start_list.append(int(subset_gff_lines["start"]))
            end_list.append(int(subset_gff_lines["end"]))

        subset_line = subset_gff_lines.iloc[0]

        element_name = str(subset_line["gene_name"])
        chrom_id = str(subset_line["seqid"])
        strandedness = str(subset_line["strand"])
        element_id = str(subset_line["id"])
        upstream = int(with_upstream)
        downstream = int(with_downstream)
        startpos = min(start_list)
        endpos = max(end_list)

        gene_occurence_dict, id_occurence_dict = extract_seq(input_fasta, chrom_dict, chrom_id, startpos, endpos,
                                                             strandedness, element_name, element_id,
                                                             upstream, downstream, gene_occurence_dict,
                                                             id_occurence_dict, output_fasta_handle, target_type,
                                                             upstream_only, downstream_only, dnaformat, with_inside,
                                                             exclude_regex)

    elif target_mode == "merged":
        for id_unique in subset_ids:
            gff_sub_lines = subset_gff_lines[subset_gff_lines["id"] == id_unique]

            startpos = list()
            endpos = list()
            for index, line in gff_sub_lines.iterrows():
                startpos.append(int(line["start"]))
                endpos.append(int(line["end"]))
                element_id = str(line["id"])
                element_name = str(line["gene_name"])
                chrom_id = str(line["seqid"])
                strandedness = str(line["strand"])
                upstream = int(with_upstream)
                downstream = int(with_downstream)

            gene_occurence_dict, id_occurence_dict = extract_seq(input_fasta, chrom_dict, chrom_id, startpos, endpos,
                                                                 strandedness, element_name, element_id,
                                                                 upstream, downstream, gene_occurence_dict,
                                                                 id_occurence_dict, output_fasta_handle, target_type,
                                                                 upstream_only, downstream_only, dnaformat, with_inside,
                                                                 exclude_regex)

    else:
        sys.exit("ERROR. Unknown target_mode.")

    return element_id, gene_occurence_dict, id_occurence_dict


def extract_seq(input_fasta, chrom_dict, chrom_id, start, end, strandedness, gene_name, gene_id, with_upstream,
                with_downstream, gene_occurence_dict, id_occurence_dict, output_fasta_handle, select_type, onlyupstream,
                onlydownstream, dnaformat, with_inside, exclude_regex, fasta_id=False):

    exclude_test = re.search(exclude_regex, gene_name)

    try:
        exclude_test.group()
        exclude_res = False
    except AttributeError:
        exclude_res = True

    if exclude_res or not exclude_regex:

        output_seq = list()

        if type(start) != list:
            start = [start]
            end = [end]

        with_upstream_orig = int(with_upstream)
        with_downstream_orig = int(with_downstream)
        with_inside_orig = int(with_inside)

        id_occurence_dict[gene_id] = id_occurence_dict.get(gene_id, 0) + 1

        if chrom_dict:
            chromosome = chrom_dict[chrom_id]

        else:
            for chrom in SeqIO.parse(input_fasta, format="fasta"):
                if chrom.id == chrom_id:
                    chromosome = chrom
                    break
            else:
                sys.exit("ERROR: Chromosome not found! Exiting...")

        chromosome_seq = chromosome.seq

        for ii in range(0, len(start)):
            element_start = int(start[ii])
            element_end = int(end[ii])

            with_inside = with_inside_orig
            if (element_end - element_start) < with_inside:
                with_inside = element_end - element_start

            if strandedness == "+":
                if ii == len(start)-1 or len(start) == 1:
                    with_downstream = with_downstream_orig
                else:
                    with_downstream = 0

                upstream = element_start - 1 - with_upstream
                if upstream < 0:
                    upstream = 0

                downstream = element_end + with_downstream
                if downstream > len(chromosome_seq):
                    downstream = len(chromosome_seq)

                if onlyupstream:
                    downstream = element_start - 1 + with_inside
                elif onlydownstream:
                    upstream = element_end

                output_seq.append(chromosome.seq[upstream:downstream])

                with_upstream = 0

            elif strandedness == "-":
                if ii == 0 or len(start) == 1:
                    with_upstream = with_upstream_orig
                else:
                    with_upstream = 0

                upstream = element_end + with_upstream

                if upstream > len(chromosome_seq):
                    upstream = len(chromosome_seq)

                downstream = element_start - 1 - with_downstream
                if downstream < 0:
                    downstream = 0

                if onlyupstream:
                    downstream = element_end - with_inside
                elif onlydownstream:
                    upstream = element_start

                output_seq_tmp = chromosome_seq[downstream:upstream]
                output_seq.append(output_seq_tmp.reverse_complement())

                with_downstream = 0

        output_seq = sum(output_seq, Seq("", generic_dna))

        fasta_description = ""
        if with_upstream_orig or with_downstream_orig:
            fasta_description += "; "
        if with_upstream_orig:
            fasta_description += str(with_upstream_orig) + " nt upstream"
        if with_upstream_orig and with_downstream_orig:
            fasta_description += " ,"
        if with_downstream_orig:
            fasta_description += str(with_downstream_orig) + " nt downstream"
        if onlyupstream:
            fasta_description += " only upstream"

            if with_inside:
                fasta_description += ", plus " + str(with_inside) + " nt inside annotation"

        if onlydownstream:
            fasta_description += " only downstream"

        gene_occurence_dict, id_occurence_dict = write_seq(output_seq, fasta_id, output_fasta_handle, dnaformat, gene_name,
                                                           gene_id, id_occurence_dict, gene_occurence_dict, select_type, fasta_description)

    return gene_occurence_dict, id_occurence_dict


def write_seq(output_seq, fasta_id, output_fasta_handle, dnaformat, gene_name, gene_id, id_occurence_dict,
              gene_occurence_dict, select_type, fasta_description):

    if fasta_id:
        gene_name = fasta_id + "." + str(id_occurence_dict[gene_id])
    elif gene_name == ".":
        gene_name = gene_id
        if select_type != "gene":
            gene_occurence_dict[gene_name] = gene_occurence_dict.get(gene_name, 0) + 1
            gene_name = gene_name + "." + str(gene_occurence_dict[gene_name])
    else:
        gene_name = gene_name + "." + gene_id + "." + str(id_occurence_dict[gene_id])

    gene_name = re.sub(" ", "_", gene_name)
    gene_name = re.sub("\t", "_", gene_name)

    if dnaformat == "single-line":
        output_fasta_handle.write(">" + gene_name + "\n")
        output_fasta_handle.write(str(output_seq) + "\n")
    else:
        SeqIO.write(SeqRecord(output_seq, id=gene_name, name="", description=fasta_description),
                    output_fasta_handle, "fasta")

    return gene_occurence_dict, id_occurence_dict


def main():
    input_fasta, input_gff, output_path, target_list, target_element, target_type, target_transcript, \
    with_upstream, with_downstream, upstream_only, downstream_only, fasta_id, \
    mychromosome, mystart, myend, mystrand, dnaformat, with_inside, exclude_raw = parse_input()

    auto_name = False

    try:
        target_transcript = int(target_transcript)
        target_transcript -= 1
        if target_transcript == -1:
            target_transcript_print = 0
        else:
            target_transcript_print = target_transcript
    except:
        target_transcript_print = target_transcript


    if output_path.lower().endswith(".fa"):
        output_fasta = output_path
    else:
        if output_path == "~/Documents/getSeq/":
            output_path = os.path.expanduser("~") + "/Documents/getSeq/"

        if output_path[-1] != "/":
            output_path += "/"

        if target_element:
            if " " in target_element:
                te_name = re.sub(" ", "-", target_element)
            else:
                te_name = target_element

            output_fasta = output_path + te_name + "." + str(target_transcript_print) + "_" + "_".join(target_type) + ".fa"

        elif target_list:
            list_name = os.path.basename(target_list)
            list_name = list_name.split(".")[0]
            output_fasta = output_path + list_name + "_" + "_".join(target_type) + ".fa"

        elif fasta_id:
            name = fasta_id.split(".")[0]
            output_fasta = output_path + name + "_userdef.fa"

        else:
            fasta_name = os.path.basename(input_fasta)
            fasta_name = fasta_name.split(".")[0]
            output_fasta = output_path + fasta_name + "_" + "_".join(target_type) + ".fa"

        auto_name = True

    try:
        dirpath = os.path.dirname(output_path)
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)
    except:
        pass

    gff_df = read_and_mod_gff(input_gff)

    if exclude_raw:
        exclude_regex = re.compile(exclude_raw)
    else:
        exclude_regex = ""

    transcript_id = find_target(input_fasta, gff_df, with_upstream, with_downstream, upstream_only, downstream_only,
                                target_list, target_element, target_type, target_transcript, output_fasta, fasta_id,
                                mychromosome, mystart, myend, mystrand, dnaformat, with_inside, exclude_regex)

    if auto_name is True and target_element:
        output_fasta_new = output_path + te_name + "." + str(transcript_id) + "_" + "_".join(target_type)

        if with_upstream:
            output_fasta_new += "_" + str(with_upstream) + "nt_up"

        if with_downstream:
            output_fasta_new += "_" + str(with_downstream) + "nt_down"

        if upstream_only:
            output_fasta_new += "_UO"
        elif downstream_only:
            output_fasta_new += "_DO"

        output_fasta_new += ".fa"

        os.rename(output_fasta, output_fasta_new)
        output_fasta = output_fasta_new

    print("Output file is " + output_fasta)

    print("Done!")


if __name__ == "__main__":
    main()

