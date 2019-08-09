import datajoint as dj
from . import get_schema_name

schema = dj.schema(get_schema_name('imaging'))


@schema
class CellType(dj.Lookup):
    definition = """
    #
    cell_type  :  varchar(100)
    ---
    cell_type_description :  varchar(4000)
    """
    contents = [
        ('Pyr', 'putative pyramidal neuron'),
        ('interneuron', 'interneuron'),
        ('PT', 'pyramidal tract neuron'),
        ('IT', 'intratelecephalic neuron'),
        ('FS', 'fast spiking'),
        ('N/A', 'unknown')
    ]


@schema
class Scan(dj.Imported):
    definition = """
    -> experiment.Session
    ---
    image_gcamp:       longblob # 512 x 512, a summary image of GCaMP at 940nm, obj-images-value1
    image_ctb:         longblob # 512 x 512, a summary image of CTB-647, obj-images-value2
    image_beads:       longblob # 512 x 512, a summary image of Beads, obj-images-value3
    recording_depth:   float    # obj-timeSeriesArrayHash-value{1, 2}-depth
    frame_time:        longblob # obj-timeSeriesArrayHash-value{1, 1}-time, aligned to the first trial start
    """

    class Roi(dj.Part):
        definition = """
        -> master
        roi_idx:   int
        ---
        -> CellType
        roi_trace:      longblob        # average fluorescence of roi obj-timeSeriesArrayHash-value{1, 1}-valueMatrix
        neuropil_trace: longblob        # average fluorescence of neuopil surounding each ROI, obj-timeSeriesArrayHash-value{1, 1}-valueMatrix
        """

@schema
class TrialTrace(dj.Imported):
    definition = """
    -> Scan.Roi
    -> experiment.SessionTrial
    ---
    trial_calcium_trace:   longblob  # trial cut using obj-timeSeriesArrayHash-value{1, 1}-trial
    """
