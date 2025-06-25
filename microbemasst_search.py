import os
import sys
import uuid

def run_microbemasst_search(
    usi: str,
    prec_tol: str = "0.05",
    mz_tol: str = "0.05",
    cos: str = "0.7",
    min_match_peaks: str = "3",
    analog_mass_below: str = "130",
    analog_mass_above: str = "120",
    use_analog: bool = False
):
    mangling = str(uuid.uuid4())
    output_temp = os.path.join("temp", "microbemasst", mangling)
    os.makedirs(output_temp, exist_ok=True)
    output_path = "../../{}".format(output_temp)
    out_file = output_path+"/fastMASST"
    print(os.getcwd())
    cmd = 'cd microbe_masst/code/ && python masst_client.py \
        --usi_or_lib_id "{}" \
        --out_file "{}" \
        --precursor_mz_tol {} \
        --mz_tol {} \
        --min_cos {} \
        --min_matched_signals {} \
        --analog_mass_below {} \
        --analog_mass_above {} \
        --database metabolomicspanrepo_index_latest \
        '.format(usi,
                 out_file,
                 prec_tol,
                 mz_tol,
                 cos,
                 min_match_peaks,
                 analog_mass_below,
                 analog_mass_above
                 )
    # Tacking on the analog flag
    if use_analog:
        cmd += " --analog true"

    print(cmd, file=sys.stderr, flush=True)
    activate_and_run = f'/bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && conda activate microbe_masst_env && {cmd}"'
    os.system("which python")
    os.system(activate_and_run)

    return output_temp

