import os
out_file_cln = open("../data/DciA_combined_cln.tsv", "a")
out_file_ipr = open("../data/DciA_combined_ipr.tsv", "a")
out_file_seqs = open("/data/scratch/janani/molevolvr_out/Dciacb_full/Dciacb.fa", "a")
with open("../data/DciA_codes.csv", "r") as infile:
    infile.readline()
    for i in range(22):
        line = infile.readline()
        line = line.rstrip().split(",")
        code = line[1]
        if line[2]:
            code = line[2]
        if line[-1] == "Yes":
            continue
        else:
            os.system("cp /data/scratch/janani/molevolvr_out/{}_full/cln_combined.tsv ../data/{}_bysp_cln_combined.tsv".format(code, code))
            fp = open("../data/{}_bysp_cln_combined.tsv".format(code), "r")
            out_file_cln.write(fp.read())
            fp.close()
            os.system("cp /data/scratch/janani/molevolvr_out/{}_full/ipr_combined.tsv ../data/{}_bysp_ipr_combined.tsv".format(code, code))
            fp = open("../data/{}_bysp_ipr_combined.tsv".format(code), "r")
            out_file_ipr.write(fp.read())
            fp.close()
            fp = open("/data/scratch/janani/molevolvr_out/{}_full/{}.fa".format(code, code))
            out_file_seqs.write(fp.read())
            fp.close()
    for line in infile:
        line = line.rstrip().split(",")
        code = line[1]
        if line[2]:
            code = line[2]
        if code == "9dxOOs":
            continue
        if line[-1] == "Yes":
            continue
        elif line[7]:
            os.system("cp /data/scratch/janani/molevolvr_out/{}_full/cln_combined.tsv ../data/{}_bydom_cln_combined.tsv".format(code, code))
            os.system("cp /data/scratch/janani/molevolvr_out/{}_full/ipr_combined.tsv ../data/{}_bydom_ipr_combined.tsv".format(code, code))
        else:
            os.system("cp /data/scratch/janani/molevolvr_out/{}_full/cln_combined.tsv ../data/{}_bysp_cln_combined.tsv".format(code, code))
            os.system("cp /data/scratch/janani/molevolvr_out/{}_full/ipr_combined.tsv ../data/{}_bysp_ipr_combined.tsv".format(code, code))

out_file_ipr.close()
out_file_cln.close()
out_file_seqs.close()