import sys, traceback, time, datetime
import pandas as pd


def decode_time(string):
    try:
        if string.count(':') == 2:
            hours, minutes, seconds = map(float, string.split(':'))
            time = datetime.timedelta(hours = hours, minutes = minutes, seconds = seconds)
        elif string.count(':') == 1:
            minutes, seconds = map(float, string.split(':'))
            time = datetime.timedelta(minutes = minutes, seconds = seconds)
        else:
            raise Exception(f"unknown time string '{string}'")
        return time.total_seconds()
    except Exception:
        traceback.print_exc()
        sys.exit("Catched exception")


def get_time_from_log(log_file_name):
    with open(log_file_name, 'r') as input_file:
        values = {"time": 0, "mem": 0}
        for line in input_file:
            if "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
                line = line.replace("Elapsed (wall clock) time (h:mm:ss or m:ss):", "").strip()
                values["time"] = decode_time(line) + values["time"]
            elif "Maximum resident set size" in line:
                values["mem"] = max(int(line.split(':')[1].strip()), values["mem"])
                print("Values", values)

        assert "time" in values, f"No time found in {log_file_name}"
        assert "mem" in values, f"No mem found in {log_file_name}"
        return values

# runtime = [
    #             get_time_from_log(
    #                 os.path.join(REPORTDIR, safe_format("safe_paths_{algorithm}", algorithm=algorithm),
    #                 "{wildcards.file_name}_k{wildcards.k}ma{wildcards.min_abundance}t{wildcards.threads}",
    #                 "log.log")
    #             )
    #             for algorithm in ALGORITHMS
    #         ]  
    
            # runtime = [get_time_from_log(os.path.join(safe_format(LOG_ALGORITHM, algorithm = algorithm))) for algorithm in FAST_ALGORITHMS],  

values = []
for log_file in snakemake.input.log_files:
    values.append(get_time_from_log(log_file))

df = pd.DataFrame(values, index=snakemake.params.row_names)
df.to_csv(snakemake.output.report, sep='\t')