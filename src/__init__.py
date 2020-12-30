from _libtimetag import *

import numpy as _np
import dateutil as _dateutil_imported

def read_sstt_header(filepath):
    """Reads the header file of a small simple time-tagged (SSTT) dataset

    Parameters
    ----------
    filepath : str
        Filepath to the header file

    Returns
    -------
    exp_header : dictionary
        A dictionary containing information describing the experiment        
    chan_header : list
        A list of dictionaries, providing information on each channel
    """
    lines = []

    header = open(filepath)

    while True:
        temp = header.readline()
        if temp == "":
            break

        temp = temp.replace('\n','')

        lines.append(temp)

    header.close()

    start_exp_header = False
    exp_header_contents = False
    start_chan_header = False
    exp_header_headings = None
    chan_header_headings = None
    chan_header_contents = False
    chan_ID_index = None

    exp_header = {'Time_unit_seconds': 81e-12,
     'device_type': 'qutau',
     'experiment_start_timestamp_UTC': None}
    chan_headers = {}

    exp_info_types = {'Time_unit_seconds': _np.double,
     'device_type': str,
     'experiment_start_timestamp_UTC': "DATETIME"}

    chan_info_types = {'ChannelID': int,
      'Filename': str,
      'NumPhotons': _np.int64,
      'NumOverflows': _np.int64,
      'Filesize': _np.int64,
      'HardwareSyncDivider': _np.int64,
      'AdditionalSyncDivider': _np.int64,
      'TotalSyncDivider': _np.int64,
      'IsPulsesChannel': _np.bool,
      'HasPulsesChannel': _np.bool,
      'CorrespondingPulsesChannel' : _np.int32,
      'HasMicrotimes' : _np.bool,
      'MicroDelayTime' : _np.int64}

    default_chan_header = {'ChannelID': None,
      'Filename': "None",
      'NumPhotons': 0,
      'NumOverflows': 0,
      'Filesize': 0,
      'HardwareSyncDivider': 1,
      'AdditionalSyncDivider': 1,
      'TotalSyncDivider': 1,
      'IsPulsesChannel': False,
      'HasPulsesChannel': False,
      'CorrespondingPulsesChannel' : None,
      'HasMicrotimes':False,
      'MicroDelayTime':0}

    for i,l in enumerate(lines):
        if l == "EXPERIMENT_HEADER":
            start_exp_header = True
            continue

        if l == "CHANNEL_HEADER":
            start_chan_header = True
            continue

        if start_exp_header:
            exp_header_headings = l.split("\t")
            start_exp_header = False
            exp_header_contents = True
            continue

        if exp_header_contents:
            contents = l.split("\t")
            if len(contents) != len(exp_header_headings):
                print("Error in experiment header!")
                return None

            for j,c in enumerate(contents):
                hd_nm = exp_header_headings[j]

                if hd_nm in exp_info_types:
                    type_ = exp_info_types[hd_nm]

                    if type_ == "DATETIME":
                        c = _dateutil_imported.parser.parse(c)
                    elif type_ == _np.double:
                        c = type_(c)
                        #c = round(c,5)
                    else:
                        c = type_(c)

                exp_header[hd_nm] = c

            exp_header_contents = False
            continue

        if start_chan_header:
            chan_header_headings = l.split("\t")

            for j,c in enumerate(chan_header_headings):
                if c == "ChannelID":
                    chan_ID_index = j
                    break


            if chan_ID_index == None:
                print("Error: could not find channel ID column")
                break

            start_chan_header = False
            chan_header_contents = True
            continue

        if chan_header_contents:
            if l == "":
                chan_header_contents = False
                continue

            contents = l.split("\t")
            contents = [c for c in contents if c]

            if len(contents) != len(chan_header_headings):
                print("Error in channel header!")
                return None

            chan_ID = int(contents[chan_ID_index])

            chan_headers[chan_ID] = default_chan_header.copy()

            for j,c in enumerate(contents):            
                # cast if possible
                header_name = chan_header_headings[j]
                if header_name in chan_info_types:
                    if chan_info_types[header_name] == _np.bool:
                        c = _np.int8(c)
                        c = _np.bool(c)

                    c = chan_info_types[header_name](c)

                    if chan_info_types[header_name] == str:
                        c = c.replace('"','')

                chan_headers[chan_ID][chan_header_headings[j]] = c
                
    return exp_header,chan_headers
            
def import_data(filepath):
    """Imports the data and header information of a small simple time-tagged (SSTT) dataset

    This function imports SSTT datasets and header information. Furthermore,
    it generates microtimes for applicable channels, i.e. photon arrival
    times relative to a reference channel, such as the laser sync channel.
    
    To only import header data, use the read_sstt_header() function.
    To only import data from one specific channel, without any preprocessing,
    use the read_sstt_data() function.
    
    Parameters
    ----------
    filepath : str
        Filepath to the header file

    Returns
    -------
    exp_header : dictionary
        A dictionary containing information describing the experiment        
    chan_header : list
        A list of dictionaries, providing information on each channel        
    data : list
        A list of dictionaries, containing the data per channel
    """
    exp_header,chan_header = read_sstt_header(filepath)
    
    data = {}

    for chan in chan_header:
        macro,micro,_ = read_sstt_data(filepath+".c"+str(chan))
        data[chan] = {}
        data[chan]["macro"] = macro
        data[chan]["micro"] = micro

    # Generate microtimes if necessary
    for chan in chan_header:            
        ch = chan_header[chan]
        
        if ch['NumPhotons'] == 0:
            ch['NumPhotons'] = len(data[chan]["macro"])
            
        if ch['HasPulsesChannel'] and not ch['HasMicrotimes'] and ch['NumPhotons'] > 0:
            # generate microtimes
            data[chan]['micro'] = gen_micro_times(data[ch["CorrespondingPulsesChannel"]]['macro'],data[chan]['macro'],chan_header[ch["CorrespondingPulsesChannel"]]["TotalSyncDivider"])
            
        if ch['IsPulsesChannel']:
            ch['PulsePeriod'] = _np.int64(round(_np.average(data[chan]['macro'][1:] - data[chan]['macro'][:-1])/ch["TotalSyncDivider"]))
        else:
            ch['PulsePeriod'] = 0
    
    return exp_header,chan_header,data