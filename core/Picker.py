import os
from datetime import timedelta as td
from pathlib import Path
from glob import glob

import seisbench.models as sbm
from obspy import read
from pandas import date_range
from obspy.core.stream import Stream

from core.PrepareData import (saveModel, loadModel, prepareWaveforms)
from core.Extra import divide_chunks
from torch import cuda
from seisbench.util import PickList


def runPicker(config):
    min_p_prob = config["min_P_probability"]
    min_s_prob = config["min_S_probability"]
    batch_size = config["batch_size"]
    overlap = config["overlap"]
    ncpu = config["number_cpu"] or os.cpu_count() - 2
    p_nam = config["picker"]
    m_nam = config["model"]
    m_upd = config["model_update"]
    m_ver = config['model_version']
    picker = eval(f"sbm.{p_nam}.from_pretrained('{m_nam}', update={m_upd}, version_str='{m_ver}')")
    print("+++ Picker loaded. Here are the picker details:",
          f"Picker: {p_nam}\nModel: {m_nam}\nupdate: {m_upd}\nversion_str: {m_ver}",
          "Doc String:",
          picker.weights_docstring,
          sep='\n')
    if cuda.is_available():
        picker.cuda()
        print("Running on GPU")
    else:
        print("Running on CPU")
    path = Path("results")
    path.mkdir(parents=True, exist_ok=True)
    startTime = config["starttime"]
    endTime = config["endtime"]
    startDateRange = date_range(startTime, endTime-td(days=1), freq="1D")
    endDateRange = date_range(startTime+td(days=1), endTime, freq="1D")
    for st, et in zip(startDateRange, endDateRange):
        print(f"+++ Run {config['picker']} Picker on period: {st} - {et}")
        dataExists = prepareWaveforms(st, et, config)
        if not dataExists:
            continue
        chunksData = glob(os.path.join("tmp", "*.mseed"))
        for c, chunkData in enumerate(divide_chunks(chunksData, 10)):
            print(f"+++ Applying SeisBench on chunk {c+1} ...")
            pick_outname = f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}_{c}"
            if not config["repick_data"] and os.path.exists(
                    os.path.join("results", f"{pick_outname}.jlib")):
                continue
            stream = Stream()
            for s in chunkData:
                stream += read(s)
            picks = picker.classify(stream,
                                    overlap=overlap,
                                    batch_size=batch_size,
                                    P_threshold=min_p_prob,
                                    S_threshold=min_s_prob,
                                    parallelism=ncpu).picks
            if config["repick_data"] and os.path.exists(
                    os.path.join("results", f"{pick_outname}.jlib")):
                os.remove(os.path.join("results", f"{pick_outname}.jlib"))

            if len(picks):
                saveModel(picks, pick_outname)
        for f in [f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.jlib",
                  f"picks_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv",
                  f"catalog_{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.csv",
                  f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}.out"]:
            if os.path.exists(os.path.join("results", f)):
                os.remove(os.path.join("results", f))
        picksList = PickList()
        picks_list = glob(os.path.join(
            "results", f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}_*.jlib"))
        picks_list = sorted(picks_list, key=lambda x: int(
            x.split("_")[-1].split(".")[0]))
        for picks in picks_list:
            picks = loadModel(picks)
            picksList += picks
        saveModel(picksList, f"{st.strftime('%Y%m%d')}_{et.strftime('%Y%m%d')}")
