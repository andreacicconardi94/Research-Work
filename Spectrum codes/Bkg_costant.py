import matplotlib.pyplot as plt
import hyperspy.api as hs
import lumispy as lum
import os, glob
import numpy as np

'''
THIS FUNCTION PERFORM A BACKGROUND CORRECTION ON THE SAMPLE
'''
def correct_bkg(s, bkg):
    #s.plot()					#Plot your sample as control
    
    ###OPTIONAL: Check the metadata of your control###
    #s.original_metadata.Object_0_channel_0.Parsed	
    
    cs = s - bkg		#Eliminate the background from the sample
    
    #cs.plot()
    
    ###WARNING: In case the minimum is lower than zero
    s = subtract_minimum(s)
    
    
    return s


'''
GRATING PHENOMENA HAPPENS BECAUSE THE SAMPLE SURFACE IS NOT PERFECTLY FLAT
WE APPLY THIS CORRECTION TO MODEL A REAL SAMPLE SURFACE WITH GRATING SURFACE
'''
def grating_correction(s):
    calibration_factor = 131072
    grating = int(s.original_metadata.Object_0_Channel_0.Parsed.SPECTROMETER.Grating__Groove_Density)

    if grating == 150:
        correction_factor_grating = 2.73E-04 # 150 gr/mm grating
    elif grating == 600:
        correction_factor_grating = 6.693659836087227e-05 # 600 gr/mm grating
    else:
        raise ImportError('Grating correction not available')
    
    fov = s.original_metadata.Object_0_Channel_0.Parsed.SITE_IMAGE.Reference_Field_of_view *1e6

    grating_calibrations = {
        'cal_factor' : calibration_factor,
        'corr_factor_grating' : correction_factor_grating,
        'field_of_view_um' : fov,
        }

    s.correct_grating_shift(*grating_calibrations.values())
    #s.plot()

    return s


def crop_n_spikes(s):
    s = crop_edges(s)
    
    #Apply this method to automathically remove spikes
    s.remove_spikes(inplace=True)
    #s.plot()
    #s.mean().plot()

    #If it does not work, you can try to remove the manually with the interactive method
    s.remove_spikes(interactive=True, inplace=True)
    s.plot()
    s.mean().plot()

    s.save("Data/No Spike.hspy")

    return s


'''
FUNCTION TO SUBTRACT A MINIMUM FROM SAMPLE AND OBTAIN A GRAPH WITH NO VALUE LOWER THAN ZERO
'''
def subtract_minimum(sample):
    minimum = int(input('Set the minimum to subtract from sample:\n'))
    l = sample - minimum
    #l.plot()
    #l.mean().plot()

    return l


'''
FUNCTION TO CROP THE EDGES OF THE MAP
'''
def crop_edges(sample):
    pixels = int(input('Insert the number of pixels to crop:\n'))
    sample.crop_edges(crop_px=pixels)

    return sample







if __name__=='__main__':
    root_folder = r'C:\Users\Cicco\Desktop\Dottorato Genova\Andrea\Sample_1_Area_02_Bkg2_3kV_100APE_1NA\Data'
    path = os.path.join(root_folder, 'HYPCard.sur')
    background_path = os.path.join(root_folder, 'BG*.txt')
    background_path = glob.glob(background_path)[0]

    sample = hs.load(path, signal_type='CL_SEM')
    background = np.loadtxt(background_path)[1]
    
    sample1_1 = correct_bkg(sample, background)
    
    sample1_2 = grating_correction(sample1_1)
    
    sample1_3 = crop_n_spikes(sample1_2)
