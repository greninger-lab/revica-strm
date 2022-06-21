import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find the accession number associated with the highest bbmap median_fold for each virus.')
    parser.add_argument('-bbmap_covstats', metavar='bbmap_covstats', required=True, type=str, help='bbmap covstats output file')
    parser.add_argument('-b', metavar='sample basename', required=True, type=str, help='sample basename')
    parser.add_argument('-m', metavar='M', type=int, default=3, help='minimum median threshold in bbmap covstats output for a reference to be considered. (Default 3)')
    parser.add_argument('-min_cov_pct', metavar='minimum covered percent for reference', type=int, default=70, help='minimum covered percent in bbmap covstats output for a reference to be considered. (Default 70)')
    args = parser.parse_args()

    init_ref_candidates = []
    inf = open(args.bbmap_covstats, 'r').readlines()
    header = inf[0]
    rec = inf[1:]

    # add to list any reference that passes the median coverage and mininum covered percent threshold
    for i in rec:
        if len(i) > 0:
            temp = i.split('\t')
            if int(temp[9]) >= args.m and float(temp[4]) >= args.min_cov_pct:
                init_ref_candidates.append(temp[0])

    # create a dictionary of references to be used for consensus calling
    init_ref_header = {}
    for i in init_ref_candidates: 
        # add .split(' ')[0] to get just the accession 
        acc = i.split(' ')[0]
        # get header tag 
        tag = i.split(' ')[1]
        # get header info
        info = " ".join(i.split(' ')[2:])
        if not tag in init_ref_header: 
            init_ref_header[tag] = [acc, info]
            
    # if reference(s) selected, output relevant info to file
    if init_ref_header:
        for i in init_ref_header:
            output_file_name = args.b + '_' + str(init_ref_header[i][0]) + '_' + i + '_' + 'vid.txt'
            output_file = open(output_file_name, 'w+')
            # output format: sample basename <tab> reference accession <tab> reference header tag <tab> reference header info
            output_file.write(args.b + "\t" + str(init_ref_header[i][0]) + "\t" + i + "\t" + str(init_ref_header[i][1]) + "\n")
            output_file.close()        

    # no reference selected, output reference with highest covered percent and read distribution info
    else:
        output_file_name = args.b + '_failed_assembly.txt'
        output_file = open(output_file_name, 'w')

        mapped_ref = []
        for line in sorted(rec, key=lambda line: float(line.split('\t')[4]), reverse=True):
            if float(line.split('\t')[4]) <=0 :
                break
            else:
                mapped_ref.append(line)

        # loop index is for: reference id, average coverage, covered percent, plus reads, minus reads, median coverage
        ref_best_cov_loop_index = [0,1,4,6,7,9]
        ref_best_cov_stats = [args.b]

        # if there are any mapped reads
        if len(mapped_ref) > 0:
            # get stats on reference with the highest covered percent
            ref_best_cov_split = mapped_ref[0].split('\t')
            for i in ref_best_cov_loop_index:
                ref_best_cov_stats.append(ref_best_cov_split[i])
            
            # calculate reads distribution based on reference tag
            mapped_stats = {}
            for i in mapped_ref:
                temp = i.split('\t')
                ref_tag = temp[0].split()[1]
                mapped_reads = int(temp[6])+int(temp[7])
                if ref_tag in mapped_stats:
                    mapped_stats[ref_tag]+=mapped_reads
                else:
                    mapped_stats[ref_tag] = mapped_reads
            mapped_stats = dict(sorted(mapped_stats.items(), key=lambda item: item[1], reverse=True))
            reads_distribution = '; '.join("{}: {}".format(k, v) for k, v in mapped_stats.items())
            ref_best_cov_stats.append(reads_distribution)

        # if no reads mapped to the database
        else:
            # +1 for reads distribution
            for i in range(0,len(ref_best_cov_loop_index)+1):
                ref_best_cov_stats.append("0")

        output_text = ""
        for i in ref_best_cov_stats:
            output_text = output_text + str(i) + "\t"

        header = "sample_name\tref_best_cov\taverage_coverage\tcovered_percent\tplus_reads\tminus_reads\tmedian_coverage\treads_distribution"
        output_file.write(header)
        output_file.write('\n')
        output_file.write(output_text)
        output_file.write('\n')
        output_file.close()
