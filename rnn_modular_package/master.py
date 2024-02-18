"""
Master script to run multiple gpu runs
use "python3 master.py -h" for help

DESCRIPTION
    
AUTHOR

    Ravi Kumar Verma

VERSION/DATE

    31st July, 2018

"""

import os
import sys
import getopt
import multiprocessing
import glob

launch_dir = os.getcwd()
sys.path.append(launch_dir)

try:
    from common_functions import add_args, write_common_log
except ImportError:
    from .common_functions import add_args, write_common_log

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

print('\nlaunch_dir =', launch_dir)


def job_submission(itr):
    """
    job submission module.
        create a command for the rnn.py
    :param itr: iterable
    """
    command = arguments
    new_kwargs = kwargs.copy()
    new_kwargs['-g'] = kwargs['-g'][itr]
    new_kwargs['-it_id'] = kwargs['-it_id'][itr]
    
    new_kwargs['-c'] = "%s_%s" % (kwargs['-c'], kwargs['-it_id'][itr])
    new_kwargs['-l'] = "%s.log" % (new_kwargs['-c'])
    
    for key, value in new_kwargs.items():
        if key != "-it_id" and value != "":
            command.append("%s %s" % (key, value))
    run_command = " ".join(str(x) for x in command)
    print('\n', run_command)
    os.system(run_command)


def show_help():
    print("\n\n\t(-:  MASTER SCRIPT  :-)\n\n")
    print("Written by:RAVI KUMAR VERMA\n")
    print("\nUSAGE:")
    print("\tpython3 master.py -r rnn_type -j jobtype -e epochs -c checkpoint_dir -i no_of_iterations -i 8 -s 0")
    print("\n\nOPTIONS:\n")
    print("\t-h or --help:\tprint help")
    print("\t-p or --pfile:\tfile with paths to input files. See path_parameters_sample.py")
    print("\t-r or --rnn:\tuse either GRU or LSTM\n\t\t\tDEFAULT: LSTM")
    print("\t-e or --epochs:\tno of optimization iterations\n\t\t\tDEFAULT: 3000")
    print("\t-c or --ckdir:\tname of checkpoint directory where checkpoints to be saved\n\t\t\tDEFAULT: checkpoints")
    print("\t-j or --jobtype: use prediction\n\t\t\tDEFAULT: prediction")
    print("\t\t\tif jobtype== validation:\tcode performs cross validation")
    print("\t\t\tif jobtype== prediction:\ttrain on whole dataset and perform analysis")
    print("\t-i or --iteration: no of iterations to be performed:")
    print("\t-s or --start_index: starting iteration")
    print("\nIMPORTANT CONSIDERATIONS:")
    print("\tUse different checkpoint directory names while doing multiple runs.")
    print("\nWORKFLOW:")
    print("\tmaster.py calls rnn.py multiple times and submit a single calculation on a given "
          "gpu\n")


def main(argv):
    if argv:
        try:
            opts, argument = getopt.getopt(argv, "hp:r:j:e:c:i:s:t:a:",
                                           ["pfile=", "rnn=", "jobtype=", "epochs=", "ckdir=", "iteration=",
                                            "start_index=", "training=", "artificial_sequences="])
        except getopt.GetoptError:
            print('\n', 'ERROR:', '\n\t', 'Check your arguments', '\n')
            show_help()
            sys.exit(2)
        
        rt = "LSTM"
        jt = "prediction"
        ep = 3000
        ck = "checkpoints"
        it = 24
        n_gpu = 8
        si = 0
        pfile = ''
        tr = ''
        afs = ''
        for opt, arg in opts:
            if opt == "-h":
                show_help()
                sys.exit()
            elif opt in ("-p", "--pfile"):
                pfile = arg
            elif opt in ("-r", "--rnn"):
                rt = arg
            elif opt in ("-j", "--jobtype"):
                jt = arg
            elif opt in ("-e", "--epochs"):
                ep = arg
            elif opt in ("-c", "--ckdir"):
                ck = arg
            elif opt in ("-i", "--iteration"):
                it = arg
            elif opt in ("-s", "--start_index"):
                si = arg
            elif opt in ("-t", "--training"):
                tr = arg
            elif opt in ("-a", "--artificial_sequences"):
                afs = arg
        
        it_ids = [int(si) + x for x in range(int(it))]
        # gpu_ids = [x for x in range(int(n_gpu))]
        gpu_ids = [1, 2, 3, 4, 5, 6, 7]
        r_kwargs = {"-g": gpu_ids, "-r ": rt, "-e": ep, "-j": jt, "-c": ck, "-it_id": it_ids, '-p': pfile, '-t': tr,
                    '-a': afs}
        
        return [pfile, r_kwargs]


if __name__ == "__main__":
    args = main(sys.argv[1:])
    parms = __import__(args[0].replace('.py', ''))
    kwargs = args[1]
    available_gpu = kwargs["-g"]
    iterations_ids = kwargs["-it_id"]
    
    clog_filename = "combined_log_itr_%s_to_%s.log" % (iterations_ids[0], iterations_ids[-1])
    arguments = ['python3', os.path.join(parms.home_dir, parms.code_dir, 'rnn.py')]
    
    """
    Multiprocessing Loop
    """
    for a in range(0, len(iterations_ids), len(available_gpu)):
        # determine how many iterations are to be run
        iterations_ids_to_run = iterations_ids[a: a + len(available_gpu)]
        iterable = [x for x in range(len(iterations_ids_to_run))]
        kwargs["-g"] = available_gpu[0: len(iterations_ids_to_run)]
        kwargs["-it_id"] = iterations_ids_to_run
        
        # ask for pool of workers
        pool = multiprocessing.Pool(processes=len(iterable))
        
        # add additional arguments obtained from commandline in the do_boot_strapping function
        func = add_args(job_submission, arguments, keywords=kwargs)
        
        # loop over iterable and submit the jobs
        pool.map(func, iterable)
        pool.close()
        pool.join()
    
    write_common_log(sorted(glob.glob("%s%s" % (kwargs["-c"], "*.log"))), clog_filename=clog_filename)
