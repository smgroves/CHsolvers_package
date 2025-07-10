import numpy as np
import pandas as pd
import itertools
import os
import re
import sys
import traceback


def read_specific_lines(file_path, line_numbers):
    result = []
    with open(file_path, "r") as file:
        for current_line_number, line in enumerate(file):
            if current_line_number in line_numbers:
                digits = []
                for i in line.strip("\n").split(","):
                    if i == "":  # at the end of each python line there is a "" appended that needs to be removed
                        pass
                    else:
                        # will do nothing if unneeded, use for python output
                        i = i.replace(" ", "")
                        i = np.real(complex(i))  # use for python output
                        digits.append(float(i))
                result.append(digits)
            if current_line_number > max(line_numbers):
                break
    result = np.array(result)
    return result

# NMG_MATLAB_2000_dt_5.50e-06_Nx_512_neumann_n_relax_4_75p_phi


def compare_L2_error_at_timepoint(indir, out_file, timepoint=None):

    # list_of_n = []
    # list_of_filenames = []

    # pattern = re.compile(fr'Nx_(\d+)_([a-zA-Z]+).*phi.csv$')
    # pattern = re.compile(fr'^([A-Za-z0-9]+)_([A-Za-z0-9]+)_(\d+).*?_Nx_(\d+).*?_(neumann|periodic).*?(?:_dtout_(\d+))?.*?_?(\d+p)?_?phi\.csv$')
    #     pattern = re.compile(
    #     r'^([A-Za-z0-9]+)_'             # Group 1: method
    #     r'([A-Za-z0-9]+)_'              # Group 2: language
    #     r'(\d+)'                        # Group 3: steps
    #     r'.*?_Nx_(\d+)'                 # Group 4: Nx
    #     r'.*?_(neumann|periodic)'       # Group 5: boundary
    #     r'(?:_dtout_(\d+))?'            # Group 6: optional dtout number
    #     r'(?:_?(25p|50p|75p))?'         # Group 7: only match these exact strings
    #     r'_phi\.csv$'                   # file ending
    # )

    # pattern = re.compile(r'^([A-Za-z0-9]+)_([A-Za-z0-9]+)_(\d+).*?_Nx_(\d+).*?_(neumann|periodic)(?:_n_relax_(4))?(?:_?(25p|50p|75p))?_phi\.csv$') #works for matlab

    # NMG_Julia_2000_dt_5.50e-06_Nx_64_periodic_dtout_1025p_phi.csv
    # should would for julia and python
    pattern = re.compile(
        rf'^([A-Za-z0-9]+)_([A-Za-z0-9]+)_(\d+).*?_Nx_(\d+).*?_(neumann|periodic)(?:_dtout_(\d+))?(?:_?(25p|50p|75p))?_?phi\.csv$')
    for root, dirs, files in os.walk(indir):
        # grab only files that match the structure
        for fname in files:
            try:
                clean_fname = re.sub(r'^.*?([A-Z]{2,4})', r'\1', fname)

                match = re.match(pattern, clean_fname)
                if match:
                    method = match.group(1)
                    language = match.group(2)
                    total_timepoints = match.group(3)
                    nx_val = match.group(4)
                    bc = match.group(5)
                    dt_out = match.group(6) or "None (10)"
                    ic = match.group(7)  # or "50p"

                    if dt_out == "None (10)":
                        dt_out_practical = 10
                    else:
                        dt_out_practical = int(dt_out)

                    # if (int(dt_out_practical) > int(timepoint)+1) & (timepoint!= 0):
                        # print(f"Timepoint not recorded. dt_out:{dt_out} is bigger than timepoint: {timepoint}")
                    # else:
                    if dt_out_practical == 2000:
                        print(
                            f"{language},{total_timepoints},{nx_val},{bc},{dt_out},{dt_out_practical},{ic} not recorded; dt_out too large.")
                        continue
                    if method == "NMG":
                        print(
                            f"{language},{total_timepoints},{nx_val},{bc},{dt_out},{dt_out_practical},{ic}")

                        NMG_fname = fname
                        SAV_fname = fname.replace("NMG", "SAV").replace("20", '2000').replace(
                            "dtout_1", "dtout_10")  # many changes for python
                        first_line = int(
                            timepoint/int(dt_out_practical) * int(nx_val))
                        last_line = int(first_line + int(nx_val))
                        line_list = range(first_line, last_line)
                        # print(f"{root}/{fname}",first_line, last_line)
                        snapshot_NMG = read_specific_lines(
                            f"{root}/{NMG_fname}", line_list)
                        snapshot_SAV = read_specific_lines(
                            f"{root}/{SAV_fname}", line_list)
                        n = np.linalg.norm(snapshot_NMG-snapshot_SAV)
                        l2_error = np.sqrt(
                            np.mean((snapshot_NMG - snapshot_SAV)**2))

                        # print(n)
                        # list_of_filenames.append(fname)
                        # list_of_n.append(n)

                        T = {}
                        T['language'] = language
                        T['GridSize'] = nx_val
                        T['boundary'] = bc
                        T['ic'] = ic
                        T['dt_out'] = dt_out
                        T["timepoint"] = timepoint
                        T['l2_error'] = l2_error
                        T['pathname'] = f"{root}/{fname}"

                        T = pd.DataFrame([T])
                        if not os.path.isfile(out_file):
                            T.to_csv(out_file, index=False)
                        else:
                            with open(out_file, "a") as f:
                                T.to_csv(f, header=False, index=False)
            except BaseException as ex:
                # Get current system exception
                ex_type, ex_value, ex_traceback = sys.exc_info()

                # Extract unformatter stack traces as tuples
                trace_back = traceback.extract_tb(ex_traceback)

                # Format stacktrace
                stack_trace = list()

                for trace in trace_back:
                    stack_trace.append("File : %s , Line : %d, Func.Name : %s, Message : %s" % (
                        trace[0], trace[1], trace[2], trace[3]))
                print(fname)
                print("Exception type : %s " % ex_type.__name__)
                print("Exception message : %s" % ex_value)
                print("Stack trace : %s" % stack_trace)
                # return list_of_filenames, list_of_n


# reset file each time if rerunning the same timepoints
out_file = f"/home/xpz5km/Cahn_Hilliard_Model_HPC/spinodal_decomp/L2_error_python.csv"


indir = "/project/g_bme-janeslab/SarahG/spinodal_decomp_06_2025/out_python"
# print("Timepoint 0")
# compare_phis_norm_at_timepoint(indir,out_file, timepoint = 0)
# print("Timepoint 20")
# compare_phis_norm_at_timepoint(indir, out_file, timepoint = 20)
# print("Timepoint 1000")
# compare_phis_norm_at_timepoint(indir, out_file, timepoint = 1000)
# print("Timepoint 2000")
# compare_phis_norm_at_timepoint(indir, out_file, timepoint = 2000)
for time in range(0, 20, 1):
    print(time)
    compare_L2_error_at_timepoint(indir, out_file, timepoint=time)
