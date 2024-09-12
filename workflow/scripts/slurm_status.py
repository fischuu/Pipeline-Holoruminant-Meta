#!/usr/bin/env python3
import re
import subprocess as sp
import shlex
import sys
import time
import logging

logger = logging.getLogger("__name__")

STATUS_ATTEMPTS = 20

jobid = sys.argv[1]

for i in range(STATUS_ATTEMPTS):
    try:
        sacct_res = sp.check_output(shlex.split("sacct -P -b -j {} -n".format(jobid)))
        res = {
            x.split("|")[0]: x.split("|")[1]
            for x in sacct_res.decode().strip().split("\n")
        }
        break
    except sp.CalledProcessError as e:
        logger.error("sacct process error")
        logger.error(e)
    except IndexError as e:
        pass
    # Try getting job with scontrol instead in case sacct is misconfigured
    try:
        sctrl_res = sp.check_output(
            shlex.split("scontrol -o show job {}".format(jobid))
        )
        m = re.search("JobState=(\w+)", sctrl_res.decode())
        res = {jobid: m.group(1)}
        break
    except sp.CalledProcessError as e:
        logger.error("scontrol process error")
        logger.error(e)
        if i >= STATUS_ATTEMPTS - 1:
            print("failed")
            exit(0)
        else:
            time.sleep(1)

status = res[jobid]

if status == "BOOT_FAIL":
    print("failed")
elif status == "OUT_OF_MEMORY":
    print("failed")
elif status.startswith("CANCELLED"):
    print("failed")
elif status == "COMPLETED":
    print("success")
elif status == "DEADLINE":
    print("failed")
elif status == "FAILED":
    print("failed")
elif status == "NODE_FAIL":
    print("failed")
elif status == "PREEMPTED":
    print("failed")
elif status == "TIMEOUT":
    print("failed")
# Unclear whether SUSPENDED should be treated as running or failed
elif status == "SUSPENDED":
    print("failed")
else:
    print("running")


# sourced from https://github.com/ManavalanG/slurm/blob/master/%7B%7Bcookiecutter.profile_name%7D%7D/slurm-status.py

# #!/usr/bin/env python3
# import subprocess
# import sys

# jobid = sys.argv[-1]

# output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

# running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED", "PREEMPTED", "TIMEOUT"]
# if "COMPLETED" in output:
#   print("success")
# elif any(r in output for r in running_status):
#   print("running")
# else:
#   print("failed")
