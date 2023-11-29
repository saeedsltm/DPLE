import warnings

from core.ExportCatalog import exporter
from core.Extra import readConfiguration
from core.FetchData import fetchRawWaveforms
from hypo71.Locate import locateHypo71
from hypocenter.Locate import locateHypocenter
from hypodd.Locate import locateHypoDD
from core.Picker import runSeisBench
from core.Visualizer import pickerStats, pickerTest, plotSeismicity
from core.CrossSection import plotCrossSection

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
        locateHypocenter(self.config)
        locateHypo71(self.config)
        locateHypoDD(self.config)

    def visualizeResults(self):
        plotSeismicity(self.config)
        pickerTest(self.config)
        pickerStats(self.config)
        plotCrossSection(self.config, "hypocenter")
        plotCrossSection(self.config, "hypo71")
        plotCrossSection(self.config, "hypoDD")


if __name__ == "__main__":
    app = Main()
    app.downloadRawData()
    app.runPicker()
    app.exportCatalog()
    app.locate()
    app.visualizeResults()
