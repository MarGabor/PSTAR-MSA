import os
import subprocess
from datetime import datetime
import traceback

def set_script_path():
    
    global SCRIPT_DIR
    SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
    os.chdir(SCRIPT_DIR)

#writing error messages to log
#void (str)
def errorFct(errorMsg):
    try:
        print(errorMsg)
        log_path = os.path.join(SCRIPT_DIR, 'job_run.log')
        with open(log_path, 'a') as errFile:
            logEntry = "\n [%s] %s" % ((datetime.now()).strftime("%d/%m/%Y %H:%M:%S"), errorMsg)
            errFile.write(logEntry)
    except:
        print(traceback.format_exc())
        open_file_err_msg = "Error log file could not be opened."
        if errorMsg == 'mp':
            return open_file_err_msg
        print(open_file_err_msg)
        exit(1)

def queue_job(job_name, output_path, job_setup_path):

    short_job_name = job_name

    if len(job_name) > 11:
        short_job_name = job_name[0:11]
        
    shell_input = []
    #shell_input = ['python3', '-m', 'PSTAR-MSA', '--wrapjob', '--msadir', os.path.join(job_setup_path, job_name), '--outputdir',
                  # output_path, '--jobname', short_job_name, '--diamondfile', '../Tools/DIAMOND/diamond', '--diamonddbfile',
                  # '../atom_diamond_db/diamond_db.dmnd', '--verbose', '--verbose', '--threads', '20', '--baseline',
                 #  '--dalibindir', '../Tools/DaliLite.v5/bin/', '--locpdbdb', '../pdb_database/', '--samplesize', '1000']

    shell_input = ['python3', '-m', 'PSTAR-MSA', '--evaljob', '--outputdir', output_path, '--jobname', short_job_name, '--verbose', '--verbose']

    #subprocess for alignment generation
    #need to PIPE stdout and stderr for output forwarding    
    process = subprocess.Popen(shell_input)

    try:
        #timeout may need to be longer (or shorter), depending on size of alignments
        #for the calculation of pairwise alignments of up to a few chains
        #500 seconds seems to be reasonable though
        #remove killing of the process, if it causes problems. Instead, add a warning.
        out, err = process.communicate()
        #if out is not None and verbosity>1:
            #print(out)
    except:
        process.kill()
        out, err = process.communicate()
        errMsg = "Communication with job subprocess failed. Job name: %s" % (job_name)
        errorFct(errMsg)
        raise

def main():

    set_script_path()

    outputs_path = '../results/'
    job_setup_path = '../jobs_to_run_no_ref/'
    
    job_name_list = []
    job_done_list = []

    #for job_done in os.listdir(outputs_path):
       # if len(job_done) > 11:
        #    job_done_list.append(job_done[0:11])
       # else:
        #    job_done_list.append(job_done)
        
    job_done_set= set(job_done_list)

    for job_name in os.listdir(job_setup_path):
        if os.path.isdir(os.path.join(job_setup_path, job_name)):
            if len(job_name) > 11:
                if job_name[0:11] in job_done_set:
                    continue
                else:
                    job_name_list.append(job_name)
            else:
                if job_name in job_done_set:
                    continue
                else:
                    job_name_list.append(job_name)

    print(job_name_list)

    for job_name in job_name_list:
        try:
            queue_job(job_name, outputs_path, job_setup_path)
        except:
            err_msg = "Exception raised during job %s. Exiting script." % job_name
            errorFct(err_msg)
            exit(1)

if __name__ == "__main__":
    main()
