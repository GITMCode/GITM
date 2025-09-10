#!/usr/bin/env python3

import io
import glob
from pathlib import Path
import os

def check_line_lengths(input_args=None):

    files = []
    glob_strs = ['*.f90', '*.f95', '*.f']

    max_line_length = 120
    warn_line_length = 100
    max_subroutine_lines = 250

    for fdir in input_args.path:
        # Inputs for `path` could be files or directories.
        # check if it's a file:
        if os.path.isfile(fdir):
            files.append(fdir)
            continue
        # Otherwise check recursively:
        for g in glob_strs:
            for path in sorted(Path(fdir).rglob(g)):
                files.append(path)

    warn_files = []
    err_files = []
    problem_lines = []
    problem_subroutine_files_and_lens = {}
    bad_subrouts = 0

    ### Read in all files, do all checks. it's not expensive.
    for f in files:

        all_lens = []
        tmp_subroutine_lengths = []
        subrout_info = []
        
        try:
            with io.open(f, 'r', encoding='utf-8') as file:
                lines = [line.rstrip() for line in file]

                for n, l in enumerate(lines):
                    all_lens.append(len(l))
                    l = l.lower()
                    if 'subroutine' in l and 'f90' in str(f):
                        end_idx = l.find('end')
                        subrout_idx = l.find('subroutine')

                        # is line commented?
                        comment_idx = l.find('!')
                        #check if the line uses ! or C for comments:
                        # if str isn't found, value is -1
                        comment_idx = comment_idx if comment_idx != -1 else  l.find('c')

                        is_commented = (comment_idx != -1) and (comment_idx < subrout_idx)
                        is_ended = (end_idx != -1) and (end_idx < subrout_idx)
                        if is_commented:
                            continue
                        try:
                            if is_ended:
                                tmp_subroutine_lengths[-1] = n - tmp_subroutine_lengths[-1]
                            elif not is_commented:
                                tmp_subroutine_lengths.append(n)
                                subrout_info.append([n, l])
                        except IndexError:
                            print(f, "\nCould not be formatted!")

        except:
            
            print(file, ' could not be opened. Ensure everything is UTF-8 encoded.')
            raise

        nLines_over_err_len = sum([1 for a in all_lens if a > max_line_length])
        nLines_over_warn_len = sum([1 for a in all_lens if a > warn_line_length])

        if nLines_over_err_len > 0:
            err_files.append(f)
            problem_lines.append(
                [l for l in range(len(all_lens)) if all_lens[l] > max_line_length])


        if nLines_over_warn_len > 0:
            warn_files.append(f)
        
        if sum([1 for a in tmp_subroutine_lengths if a > max_subroutine_lines]) > 0:
            problem_subroutine_files_and_lens[str(f)] = []
            for n, ilen in enumerate(tmp_subroutine_lengths):
                if ilen > max_subroutine_lines:
                    problem_subroutine_files_and_lens[str(f)].append([subrout_info[n], ilen])
                    bad_subrouts += 1
        
    # Always show errors (lines longer than max_line_length characters):
    if len(err_files) > 0:
        print("\n\n=============\n\n")
        for n,(f,l) in enumerate(zip(err_files, problem_lines)):
            print(f"{f}\n\tLines: {l}")
        print("\n\n=============\n\n")
        

        raise ValueError(
            f"""\n >>>> ERROR!! \n    >> LINES TOO LONG!
    >> The {len(err_files)} file(s) above contain lines longer than {max_line_length} characters.
    >> These must be fixed before your code will be accepted. Please fix and try again.\n""")


    # only show warnings (above recommended) if the users asks.
    if len(warn_files) > 0 and input_args.line_length:

        print(f"""
>> WARNING: {len(warn_files)} files have lines longer than the maximum 
    recommended line length of {warn_line_length} characters. 
    It is recommended, but not necessary to fix this.""")


    # Only show subroutine length if the user asks:
    if bad_subrouts > 0 and input_args.subroutines:
        if input_args.verbose:
            for fname, info_len in problem_subroutine_files_and_lens.items():
                print(f"\n\n  within file: {fname}:")
                for i in info_len:
                    print(f"line: {i[0]} has length {i[-1]}")

        print(f"""
>> WARNING: {bad_subrouts} subroutines are over the recommended length of
    {max_subroutine_lines} lines. 
    It is recommended, but not necessary to refactor these.
      """)



def check_git_status(args):
    import subprocess

    gitstatus = subprocess.run('git status --porcelain', shell=True, capture_output=True)
    if len(gitstatus.stdout.decode()) == 0:
        print('git status appears OK! ')
        return
    
    if args.verbose:
        gitstatus = gitstatus.stdout.decode().split('\n')
        print("Git status not passed for files:\n", gitstatus)
    print('done')
    return

def find_fpretty_config_file(args):
    """
    Fprettify doesn't run right without the config file. This will try to find one.

    check:
        - input fpretty_config argument
        - $pwd
        - $home
        - target directory
        - target_directory/srcPython/fprettify/
    
    """
    fpretty_opts = []

    if args.fpretty_config:
        if os.path.isfile(args.fpretty_config):
            return args.fpretty_config
        else:
            print("==> Invalid config file. Looking elsewhere...")
    fpretty_opts.append(os.getcwd())
    fpretty_opts.append(os.path.expanduser("~"))
    fpretty_opts.extend(args.path)
    fpretty_opts.extend([os.path.join(p, 'srcPython/fprettify') for p in args.path])

    for test_dir in fpretty_opts:
        test_fpretty_file = os.path.join(test_dir, '.fprettify.rc')
        if os.path.isfile(test_fpretty_file):
            print(f"===> Using fprettify config file '{test_fpretty_file}' ")
            print("     Specify --fpretty_config [path] if you wish to use a different file.")
            return test_fpretty_file
        
    raise FileNotFoundError("Fprettify config file not found in $HOME, pwd, or target directory!")


def run_fprettify(args):
    import subprocess
    try: #check if fprettify is installed.
        import fprettify
    except ModuleNotFoundError:
        print("==> fprettify installation not found! <==\n")
        print("\tGo to https://github.com/GITMCode/fprettify to install.\n")
        raise

    # Find the absolute path to the config file, set the string so its easier to call later
    config_file_path = find_fpretty_config_file(args)

    config_str = " -c " + config_file_path

    if args.format_check:
        print('Checking format with fprettify...')

        fpretty_out = subprocess.run(
            f"fprettify {config_str} -d {' '.join(args.path)}",
            shell=True, capture_output=True)

        if args.verbose:
            print(fpretty_out.stdout.decode())

        if len(fpretty_out.stdout) > 0:
            raise ValueError("Format does not comply."
            "\n==> Format check NOT passed. Run with -v/--verbose to see proposed changes."
            "\n    Exiting. Please verify all files match standards and retry.")

    elif args.auto_format:
        print('Auto-formatting with fprettify... Changes will not be listed.')

        fpretty_out = subprocess.run(
            f"fprettify {config_str} {' '.join(args.path)}",
            shell=True, capture_output=True)
    



def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="")

    parser.add_argument(
        "path",
        type=str,
        default=['.'],
        nargs="*",
        help="Input file(s) or directories.\n"
             "If the input is a directory all files with extension "
             "f90 and f95 are checked."
    )

    parser.add_argument(
        "-g", "--git_status", action="store_true", default=False,
        help="Show output of git status? (Default: False)"
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Show detailed information on which files cause warnings?"
    )

    parser.add_argument(
        "-s", "--subroutines", action="store_true", default=False,
        help="Show warnings for subroutines which are too long? (Default: False)"
    )

    parser.add_argument(
        "-l", "--line_length", action="store_true", default=False,
        help="Show warnings for lines above recommended length? (Default: False)"
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-a", "--auto_format", action="store_true", default=False,
        help="Attempt to automatically format Fortran files? (Default: False)"
        " This will make changes in-place."
    )

    group.add_argument(
        "-f", "--format_check", action="store_true", default=False,
        help="Check if input files are formatted correctly?"
        " Errors will be raised if `fprettify` recommends any changes. (Default: False)"
    )

    parser.add_argument(
        "--fpretty_config", default=".fprettify.rc", type=str,
        help="Path to fprettify configuration file. (default: '.fprettify.rc')"
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    print("""=================================================================
This script will aid you in formatting your code before 
submitting a pull request to https://github.com/GITMCode/GITM.git

It is by no means exhaustive however, and you should double 
check all changes made. Use -h/--help flags to see options.
================================================================
""")

    input_args = parse_arguments()

    if input_args.git_status:
        check_git_status(input_args)
    

    if input_args.auto_format or input_args.format_check:
        check_line_lengths(input_args)
        run_fprettify(input_args)
        # check line lengths again in case fprettify messed anything up
        try:
            check_line_lengths(input_args) 
        except:
            print(f"==> Original line length test and fprettify appear to have run successfully.\n"
                  "   Fprettify auto-formatting appears to have made lines too long.\n"
                  "   Manual inspection of the changed files is required\n"
            )
            raise #this should raise the correct errors... we'll see!
    else:
        check_line_lengths(input_args)


    print('No critical warnings generated.')
