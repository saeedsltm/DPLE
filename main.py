from core.Extra import readConfiguration
from core.FetchData import fetchRawWaveforms
from core.Picker import runSeisBench
from core.Visualizer import plotSeismicity, pickerTest, pickerStats
from core.ExportCatalog import exporter
from hypocenter.Locate import locateHypocenter
from hypo71.Locate import locateHypo71
from hypodd.Locate import locateHypoDD
import warnings
warnings.filterwarnings("ignore")


class Main():

    def __init__(self):
        self.config = readConfiguration()

    def downloadRawData(self):
        fetchRawWaveforms(self.config)

    def runPicker(self):
        runSeisBench(self.config)

    def exportCatalog(self):
        exporter(self.config)

    def locate(self):
        # locateHypocenter(self.config)
        # locateHypo71(self.config)
        locateHypoDD(self.config)

    def visualizeResults(self):
        plotSeismicity(self.config)
        pickerTest(self.config)
        pickerStats(self.config)


if __name__ == "__main__":
    app = Main()
    # app.downloadRawData()
    # app.runPicker()
    # app.exportCatalog()
    app.locate()
    # app.visualizeResults()

