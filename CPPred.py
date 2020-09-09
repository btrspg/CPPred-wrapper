#! /usr/bin/env python3

import argparse as agp
import os
import tempfile
from subprocess import check_call

import Bio.SeqIO as Seq

from cppred import CTD
from cppred import FrameKmer
from cppred import ORF_length as len
from cppred import ProtParam as PP
from cppred import fickett
from cppred.paras import get_model_range_hexamer


def get_tempdir():
    return tempfile.TemporaryDirectory()


def coding_nocoding_potential(input_file):
    coding = {}
    noncoding = {}
    for line in open(input_file).readlines():
        fields = line.split()
        if fields[0] == 'hexamer': continue
        coding[fields[0]] = float(fields[1])
        noncoding[fields[0]] = float(fields[2])
    return coding, noncoding


def output_feature(seq_file, hex_file, species, tmpdir):
    tmp = open(path_file(tmpdir, 'test.f_svm'), 'w')
    feature = open(path_file(tmpdir, 'test.feature'), 'w')
    out_label = 1
    coding, noncoding = coding_nocoding_potential(hex_file)
    if species == "Human":
        feature.write("\t".join(map(str, ["#ID", "ORF-integrity", "ORF-coverage", "Instability", "T2", "C0", "PI",
                                          "ORF-length", "AC", "T0", "G0", "C2", "A4", "G2", "TG", "A0", "TC", "G1",
                                          "C3", "T3", "A1", "GC", "T1", "G4", "C1", "G3", "A3", "Gravy", "Hexamer",
                                          "C4", "AG", "Fickett", "A2", "T4", "C", "G", "A", "T"])) + "\n")
    if species == "Integrated":
        feature.write("\t".join(map(str,
                                    ["#ID", "ORF-coverage", "ORF-integrity", "GC", "Instability", "ORF-length", "T0",
                                     "Fickett", "G2", "C3", "PI", "A3", "C1", "G3", "Hexamer", "TG", "G1", "TC", "A0",
                                     "A1", "AC", "C2", "G0", "T4", "C0", "A4", "G", "A2", "T", "T3", "G4", "C4",
                                     "Grary", "T2", "AG", "AT", "T1", "A", "C"])) + "\n")
    for seq in Seq.parse(seq_file, 'fasta'):
        seqid = seq.id
        A, T, G, C, AT, AG, AC, TG, TC, GC, A0, A1, A2, A3, A4, T0, T1, T2, T3, T4, G0, G1, G2, G3, G4, C0, C1, C2, C3, C4 = CTD.CTD(
            seq.seq)
        insta_fe, PI_fe, gra_fe = PP.param(seq.seq)
        fickett_fe = fickett.fickett_value(seq.seq)
        hexamer = FrameKmer.kmer_ratio(seq.seq, 6, 3, coding, noncoding)
        Len, Cov, inte_fe = len.len_cov(seq.seq)
        if species == "Human":
            tem = [inte_fe, Cov, insta_fe, T2, C0, PI_fe, Len, AC, T0, G0, C2, A4, G2, TG, A0, TC, G1, C3, T3, A1, GC,
                   T1, G4, C1, G3, A3, gra_fe, hexamer, C4, AG, fickett_fe, A2, T4, C, G, A, T]
            feature.write("\t".join(map(str,
                                        [seqid, inte_fe, Cov, insta_fe, T2, C0, PI_fe, Len, AC, T0, G0, C2, A4, G2, TG,
                                         A0, TC, G1, C3, T3, A1, GC, T1, G4, C1, G3, A3, gra_fe, hexamer, C4, AG,
                                         fickett_fe, A2, T4, C, G, A, T])) + "\n")
        elif species == "Integrated":
            tem = [Cov, inte_fe, GC, insta_fe, Len, T0, fickett_fe, G2, C3, PI_fe, A3, C1, G3, hexamer, TG, G1, TC, A0,
                   A1, AC, C2, G0, T4, C0, A4, G, A2, T, T3, G4, C4, gra_fe, T2, AG, AT, T1, A, C]
            feature.write("\t".join(map(str,
                                        [seqid, Cov, inte_fe, GC, insta_fe, Len, T0, fickett_fe, G2, C3, PI_fe, A3, C1,
                                         G3, hexamer, TG, G1, TC, A0, A1, AC, C2, G0, T4, C0, A4, G, A2, T, T3, G4, C4,
                                         gra_fe, T2, AG, AT, T1, A, C])) + "\n")
        else:
            raise ValueError('Species could only be Human or Integrated, not ' + species)
        tmp.write(str(out_label) + '\n')
        for label, item in enumerate(tem):
            tmp.write(str(label + 1) + ':' + str(item) + '\n')
        tmp.write('\n')
    tmp.close()


def print_and_run(cmd: str):
    print("CMD:" + cmd)
    check_call(cmd.split(), shell=True)


def path_file(path, f):
    return os.path.join(path, f)


def predict(range_file, model_file, libsvm_bin, tmpdir):
    svm_scale = path_file(libsvm_bin, 'svm-scale') + ' -r ' + range_file + ' ' + \
                ' -s ' + path_file(tmpdir, 'test.scaled ') + path_file(tmpdir, 'test.f_svm ')

    # os.system(libsvm_bin + '/svm-scale -r ' + range_file + ' test.f_svm  > test.scaled ')
    # os.system(svm_scale)
    print_and_run(svm_scale)
    svm_preict = path_file(libsvm_bin, 'svm-predict ') + ' -b 1 ' + path_file(tmpdir, 'test.scaled') + \
                 ' -s ' + path_file(tmpdir, 'tmp2.txt ') \
                 + model_file + path_file(tmpdir, 'tmp.txt ')
    # os.system(libsvm_bin + '/svm-predict -b 1 test.scaled ' + model_file + ' tmp.txt >  tmp2.txt')
    print_and_run(svm_preict)

    coding_poten = open(path_file(tmpdir, 'coding_potential'), 'w')
    coding_poten.write("\t".join(map(str, ["table", "coding_potential"])) + "\n")

    for line in open(path_file(tmpdir, 'tmp.txt'), 'r').readlines():
        if line[0] == "l":
            continue
        coding_potential = line.split(" ")[1]
        if line.split(" ")[0] == "1":
            coding_poten.write("\t".join(map(str, ["coding", coding_potential])) + "\n")
        else:
            coding_poten.write("\t".join(map(str, ["noncoding", coding_potential])) + "\n")


def merge(tmpdir, output_file):
    test_feature = path_file(tmpdir, 'test.feature')
    coding_potential = path_file(tmpdir, 'coding_potential')
    cmd = 'paste {tf} {cp} > {outfile}'.format(
        tf=test_feature,
        cp=coding_potential,
        outfile=output_file
    )
    print_and_run(cmd)
    # print("CMD:" + "paste test.feature coding_potential >" + output_file)
    # os.system("paste test.feature coding_potential >" + output_file)


def deleted(tmpdir):
    tmpdir.cleanup()
    # os.system("rm test.*")
    # os.system("rm coding_potential")
    # os.system("rm tmp*")


def main():
    parser = agp.ArgumentParser()
    parser.add_argument('-i', '--RNA_file', help="the input FASTA file of RNA sequence")
    parser.add_argument('-hex', '--hexmar', help="the input of hexmar", default=None)
    parser.add_argument('-r', '--range', help="the input file of range", default=None)
    parser.add_argument('-mol', '--model', help="the input file of model", default=None)
    parser.add_argument('-spe', '--species', help="the input species", required=True,
                        choices=['Human', 'Integrated'], default='Human')
    parser.add_argument('-o', '--outfile', help="output file")
    parser.add_argument('--libsvm-bin', default='/usr/local/bin', help='libsvm bin')
    args = parser.parse_args()
    tmp = get_tempdir()
    print("temporary direction: " + tmp.name)
    m_f, r_f, h_f = get_model_range_hexamer(args.species)
    if None is not args.hexmar:
        h_f = args.hexmar
    if None is not args.model:
        m_f = args.model
    if None is not args.range:
        r_f = args.range
    output_feature(args.RNA_file, h_f, args.species, tmp.name)
    os.system('head ' + path_file(tmp.name, 'test.f_svm'))
    predict(r_f, m_f, args.libsvm_bin, tmp.name)
    merge(tmp.name, args.outfile)
    deleted(tmp)


if __name__ == '__main__':
    main()
