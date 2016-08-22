import math
import numpy as np
import numpy, scipy.sparse
from operator import itemgetter
import glob,re,argparse,datetime,os,os.path,json
import shlex
from subprocess import Popen, PIPE
from multiprocessing import Pool

def my_print(s):
    print("[" + str(datetime.datetime.now()) + "] " + s)

parser = argparse.ArgumentParser(description='Experiments for "FIB Efficiency in Distributed Platforms" (ICNP 2016).')
parser.add_argument('-d', '--data', default='data', help="Input FIB folder; the FIBs should be there in .txt format.")
parser.add_argument('-o', '--out', default='out.tex', help="Output file where to print the tables in LaTeX format.")
parser.add_argument('-t', '--threads', default=4, help="Number of threads to use for FIB processing.")
args = parser.parse_args()

def divide_fib(fname_in, fname_out_tx, fname_out_rx_0, fname_out_rx_1):
    acts = {}
    rows = []
    with open(fname_in) as f:
        for line in f:
            arr = line.strip().split()
            acts[arr[1]] = acts.get(arr[1], 0) + 1
            rows.append(arr)
    acts_sorted = sorted([(k,v) for k,v in acts.items()], key=itemgetter(1), reverse=True)
    if acts_sorted[0][1] > (len(rows) / 2.5):
        acts_first = { acts_sorted[1][0] : True }
        acts_second = { k[0] : True for k in acts_sorted[2:]}
        num_from_firstaction = (len(rows) / 2) - acts_sorted[1][1]
    else:
        i = 0
        for i in range(len(acts_sorted)):
            if np.sum([x[1] for x in acts_sorted[:i]]) > len(rows) / 2:
                break
        acts_first = { k[0] : True for k in acts_sorted[:i]}
        acts_second = { k[0] : True for k in acts_sorted[i:]}
        num_from_firstaction = acts_sorted[0][1]

    cur_from_firstaction = 0
    with open(fname_out_tx, "w") as fout_divided:
        with open(fname_out_rx_0, "w") as fout_egress_first:
            with open(fname_out_rx_1, "w") as fout_egress_second:
                for row in rows:
                    if row[1] == acts_sorted[0][0]:
                        if cur_from_firstaction > num_from_firstaction:
                            res = 1
                        else:
                            res = 0
                        cur_from_firstaction += 1
                    elif row[1] in acts_first:
                        res = 0
                    else:
                        res = 1
                    fout_divided.write("%s\t%d\n" % (row[0], res) )
                    if res == 0:
                        fout_egress_first.write("%s\t%s\n" % (row[0], row[1]) )
                    else:
                        fout_egress_second.write("%s\t%s\n" % (row[0], row[1]) )

def get_exitcode_stdout_stderr(cmd):
    """
    Execute the external command and get its exitcode, stdout and stderr.
    """
    args = shlex.split(cmd)
    proc = Popen(args, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    exitcode = proc.returncode
    return exitcode, out, err

def print_one_table(outf, res_ks, comment, value_func):
    outf.write("%% %s\n\n" % comment)
    for k in res_ks:
        outf.write(" % ".join(value_func(results[k])) + " \\\\\n")
    outf.write("\n\n")

def get_omr(v, prefix, egress=True):
    res = [
        v[prefix + " 1gr wid"],
        v[prefix + " 1gr I"],
        v[prefix + " 1gr D"],
        v[prefix + " 1gr size"],
        v[prefix + " #gr"],
        v[prefix + " size"]
    ]
    if egress:
        res += [ v[prefix + " egress size"] ]
    return res

def process_one_file(in_fname):
    my_print("%s" % in_fname)
    if (not os.path.exists(in_fname[:-4] + '.TX.div.tmp')) or (not os.path.exists(in_fname[:-4] + '.RX.0.tmp')) or (not os.path.exists(in_fname[:-4] + '.RX.1.tmp')):
        my_print("\t...divided egress FIBs do not (all) exist, recreating...")
        divide_fib(in_fname, in_fname[:-4] + '.TX.div.tmp', in_fname[:-4] + '.RX.0.tmp', in_fname[:-4] + '.RX.1.tmp')
    if not os.path.exists(in_fname[:-4] + '.res.json'):
        my_print("\t...%s.res.json does not exist, running experiments..." % in_fname[:-4])
        cmd = "./bin/reduce_fib -i %s -r %s.TX.div.tmp -0 %s.RX.0.tmp -1 %s.RX.1.tmp" % (in_fname, in_fname[:-4], in_fname[:-4], in_fname[:-4])
        exitcode, out, err = get_exitcode_stdout_stderr(cmd)
        result = {}
        for line in out.decode('utf8').split('\n'):
            arr = line.strip().split('\t')
            if len(arr) == 2:
                result[arr[0]] = arr[1]
        with open("%s.res.json" % in_fname[:-4], 'w') as outf:
            json.dump(result, outf)
        my_print("\t...experiments done and dumped to %s.res.json." % in_fname[:-4])

if __name__ == "__main__":
    out_file = args.out
    threadpool = Pool(args.threads)


    input_dir = "data_test" # args.data
    my_print("Processing FIBs from directory %s/ in %d threads..." % (input_dir, args.threads) )
    threadpool.map(process_one_file, glob.glob('%s/*.txt' % input_dir))

    my_print("All FIBs processed, reading results and compiling LaTeX tables in %s..." % out_file)
    results = {}
    for in_fname in glob.glob('%s/*.txt' % input_dir):
        results[in_fname[:-4].split('/')[-1]] = json.load( open("%s.res.json" % in_fname[:-4]) )
    res_ks = sorted(list(results.keys()))
    with open(out_file, 'w') as outf:
        print_one_table(outf, res_ks, "Table I. One-stage forwarding (Fig. 2a) and two-stage forwarding (Fig. 2b).",
            lambda v : [
                v["Original actions"], v["Original rules"], v["Original size"],
                v["1-stage Boolean rules"], v["1-stage Boolean size"],
                v["2-stage ingress orig rules"], v["2-stage ingress orig size"],
                v["Egress 0 Bool rules"] + v["Egress 1 Bool rules"],
                v["Egress 0 Bool size"] + v["Egress 1 Bool size"]
                ])

        print_one_table(outf, res_ks, "Table II. Equivalent FIB representations: one-stage ingress, filter-order ind., false positive check on egress (Fig. 3a); two-stage forwarding, RX action-order ind., TX lookup (Fig. 3b); two stage forwarding, RX non-confl. rules and TX lookup (Fig. 3c).",
            lambda v :
                get_omr(v, "1-stage RX filter, TX fpcheck") + 
                get_omr(v, "2-stage RX action") + 
                get_omr(v, "2-stage RX nonconfl")
            )

        print_one_table(outf, res_ks, "Table III. Non-equivalent FIB representations: one-stage, RX action-order ind., no check on egress; (Fig. 4a); two-stage, RX action-order ind., TX action-order ind. lookup; (Fig. 4b); two-stage, RX non-confl. rules, TX action-order ind. lookup (Fig. 4c).",
            lambda v : 
                get_omr(v, "1-stage RX action, TX none") + 
                get_omr(v, "2-stage RX action", egress=False) + 
                [ str(float(v["Divided egress 0 size"]) + float(v["Divided egress 1 size"])) ] +
                get_omr(v, "2-stage RX nonconfl", egress=False) +
                [ str(float(v["Divided egress 0 size"]) + float(v["Divided egress 1 size"])) ]
            )

        print_one_table(outf, res_ks, "Table IV. Summary table: ingress, egress, and total classifier size for all FIB representations.",
            lambda v :  [
                v["1-stage Boolean size"],
                v["2-stage ingress orig size"],
                v["Total divided egress size"],
                str(float(v["2-stage ingress orig size"]) + float(v["Total divided egress size"])),
                v["1-stage RX filter, TX fpcheck size"],
                v["1-stage RX filter, TX fpcheck egress size"],
                str(float(v["1-stage RX filter, TX fpcheck size"]) + float(v["1-stage RX filter, TX fpcheck egress size"])),
                v["2-stage RX action size"],
                v["2-stage RX action egress size"],
                str(float(v["2-stage RX action size"]) + float(v["2-stage RX action egress size"])),
                v["2-stage RX nonconfl size"],
                v["2-stage RX nonconfl egress size"],
                str(float(v["2-stage RX nonconfl size"]) + float(v["2-stage RX nonconfl egress size"])),
                str(float(v["1-stage RX action, TX none size"]) + float(v["1-stage RX action, TX none egress size"])),
                v["2-stage RX action size"],
                str(float(v["Divided egress 0 size"]) + float(v["Divided egress 1 size"])),
                str(float(v["2-stage RX action size"]) + float(v["Divided egress 0 size"]) + float(v["Divided egress 1 size"])),
                v["2-stage RX nonconfl size"],
                str(float(v["Divided egress 0 size"]) + float(v["Divided egress 1 size"])),
                str(float(v["2-stage RX nonconfl size"]) + float(v["Divided egress 0 size"]) + float(v["Divided egress 1 size"])),
            ])


