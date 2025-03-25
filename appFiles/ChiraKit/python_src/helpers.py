
import numpy  as np
import pandas as pd
import re

temperature_to_kelvin = lambda T: T + 273.15 if np.max(T) < 270 else T

def solve_one_root_quadratic(a,b,c):

    return 2*c / (-b - np.sqrt(b**2 - 4*a*c))

def solve_one_root_depressed_cubic(p,q):
    # X**3 + pq + q = 0

    delta = np.sqrt((q**2/4) + (p**3/27))

    return np.cbrt(-q/2+delta) + np.cbrt(-q/2-delta)

def concat_signal_lst(listOfSignals):
    
    try:
        allSignal       = np.concatenate(listOfSignals)
    except:
        allSignal       = listOfSignals

    return(allSignal)

def round_to_significant_digits(arr, digits):
    rounded_arr = np.zeros_like(arr)
    for i, value in np.ndenumerate(arr):
        if np.abs(value) < np.finfo(float).eps:
            rounded_arr[i] = value  # Keep small numbers close to zero as is
        else:
            magnitude = np.power(10, digits - np.floor(np.log10(np.abs(value))) - 1)
            rounded_arr[i] = np.around(value * magnitude) / magnitude
    return rounded_arr

def string_to_units(string):

    '''
    Convert the string into one of the followings:

        'milliabsorbance', 'molarExtinction', 'degrees', 'millidegrees', 'molarEllipticity',
        'meanUnitMolarEllipticity', 'meanUnitMolarExtinction'

    Useful function to reverse the output from the function 'workingUnits2ProperLabel' from 'helpers_plotting.R'

    E.g. 
    Input - 'Millidegrees (m°)'      /   Output - millidegrees
    Input - 'Milliabsorbance (mΔA)'  /   Output - milliabsorbance
    Input - 'Δε ...'                /   Output - meanUnitMolarExtinction

    '''

     # Set milidegrees as default output
    units = 'millidegrees'

    string = string.lower()

    unit_mappings = {
    'absorbance'            : 'absorbance',
    'degrees'               : 'degrees',
    'mue'                   : 'meanUnitMolarEllipticity',
    'δε'                    : 'meanUnitMolarExtinction',
    'delta epsilon'         : 'meanUnitMolarExtinction',
    'molar extinction'      : 'molarExtinction',
    'molar ellipticity'     : 'molarEllipticity'
    }

    for key, value in unit_mappings.items():
        if key in string:
            units = value
            break  # Exit the loop once a match is found

    if 'milli' in string or 'mili' in string:

        units = 'milli' + units

    return units

def guess_input_units_from_metadata_dictionary(metadata_dic):

    for key, value in metadata_dic.items():

        key_lower = key.lower()

        # Known case, exported files from ChiraKit
        if 'Units of the CD signal' in key:

            return string_to_units(value)

    # Set 'Millidegrees' as default output
    return 'millidegrees'


def guess_parameter_from_metadata_dictionary(metadata_dic,parameter_names):

    '''
    parameter_names should be a list of strings 
    E.g., ['concentration','conc'], 
    '''

    for key, value in metadata_dic.items():

        key_lower = key.lower()

        for parameter in parameter_names:

            if parameter in key_lower:

                try:

                    return float(value)

                except:

                    pass

    return 0

def filter_matrix_by_vector(np_matrix,np_vector,min_value,max_value):

    """
    Filter the matrix using a certain range

    Requires: 
        
        1) The matrix 'np_matrix' of dimensions
    n*m where n matches the length of np_vector 
        2) The vector 'np_vector' of length n
        3) The lower bound 'min_value'
        4) The upper bound 'max_value'

    Returns the filtered matrix 
    """

    np_tog = np.column_stack((np_matrix, np_vector))
    tot    = np_matrix.shape[1]
    np_tog = np_tog[np_tog[:,tot]  >= min_value]
    np_tog = np_tog[np_tog[:,tot]  <= max_value]
    np_tog = np.delete(np_tog, tot, 1)

    return(np_tog)

def filter_vector_by_values(np_vector,min_value,max_value):

    """
    Filter the vector using a certain range

    Requires: 
        
        1) The vector 'np_vector' of length n.
        2) The lower bound 'min_value'
        3) The upper bound 'max_value'


    Returns the filtered vector 
    """

    np_temp = np_vector[np_vector   >= min_value]
    np_temp = np_temp[np_temp       <= max_value]

    return(np_temp)

def get_temp_at_maximum_of_derivative(temps,signal_derivative):

    tms_derivative = np.take(temps, np.argmax(signal_derivative,axis=0)) 
    return tms_derivative

def get_temp_at_minimum_of_derivative(temps,signal_derivative):

    tms_derivative = np.take(temps, np.argmin(signal_derivative,axis=0)) 
    return tms_derivative

def extract_words(input_string):

    # Pattern to match words with only alphabet characters or underscores (without hyphens)
    pattern = r'\b[A-Za-z_]+\b'

    words = re.findall(pattern, input_string)

    return words

def clean_function_text(text):

    # Use numpy notation
    cleaned_text = text.replace('e^', 'exp')
    cleaned_text = cleaned_text.replace('exp(','np.exp(') 
    cleaned_text = cleaned_text.replace('^', '**') 
    cleaned_text = cleaned_text.replace('log(','np.log(')
    cleaned_text = cleaned_text.replace('sqrt(','np.sqrt(')

    return cleaned_text

def signal_at_222nm(signal,wavelength):

    '''
    Extract the signal at 222 nm from the CD spectra
    '''

    if 222 in wavelength:

        signal_222nm = signal[wavelength == 222,:]
        return signal_222nm.flatten()

    signal_222nm = []

    signal_temp     = signal[wavelength >= 221,:]
    wavelength_temp = wavelength[wavelength >= 221]

    signal_temp     = signal_temp[wavelength_temp <= 223,:]
    wavelength_temp = wavelength_temp[wavelength_temp <= 223]

    # Sort in increasing order
    # Get the indices that would sort the wavelength vector
    sorted_indices = np.argsort(wavelength_temp)

    # Use the sorted_indices to rearrange the rows of the matrix and vector
    signal_temp     = signal_temp[sorted_indices,:]
    wavelength_temp = wavelength_temp[sorted_indices]

    # Total CD spectra of this experiment
    n_spectra = signal.shape[1]

    # Iterate over the columns
    for j in range(n_spectra):
        # Interpolate
        y = np.interp(222, wavelength_temp, signal_temp[:, j], left=None, right=None, period=None)

        signal_222nm.append(y)

    return np.array(signal_222nm)

## check or remove code below!!!
def check_good_parameters_estimation(params,low_bound,high_bound,params_name):

    """
        
    Check that the estimated parameters are far from the bounds, 
    if not there is a problem with the fitting.

    For the first  4 params 'kN', 'bN', 'kU', 'bU' we will normalize low_bound - high_bound range to 0-1
    and then verify that the parameters lie in the interval 0.02-0.98.
        
    """

    low_bound  = np.array(low_bound)
    high_bound = np.array(high_bound)
    params     = np.array(params)

    params_normalized = (params[:4] - low_bound[:4] ) / (high_bound[:4] - low_bound[:4])

    lie_in_correct_interval   = np.logical_and(params_normalized < 0.98,
        params_normalized > 0.02)

    """

    For Tm or T1 and T2 and T_onset or T_onset1 and T_onset2  we will check that they are 1 degree from the boundaries

    """

    for temp_param_name in ["Tm","T1","T2","T_onset","T_onset1","T_onset2"]:

        if temp_param_name in params_name:

            position_of_tm = params_name.index(temp_param_name) 
            tm             = params[position_of_tm]
            tm_bound_low   = low_bound[position_of_tm]
            tm_bound_up    = high_bound[position_of_tm]

            lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([tm > (tm_bound_low+1) and tm < (tm_bound_up-1)]))

    """

    For Tm or dHm, dHm1 and dHm2 we will check that they are between 2.5 and 750 kcal/mol. 10466 and 3098230 in Joules

    """

    for dh_param_name in ["dHm","dHm1","dHm2"]:

        if dh_param_name in params_name:

            position_of_dh = params_name.index(dh_param_name) 
            dh             = params[position_of_dh]
            lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([dh > 10466 and dh < 3098230]))

    """

    For Ki we will check that it lies between 1e-3 and 1e3

    """

    if "Ki" in params_name:
        position_of_ki = params_name.index("Ki") 
        ki             = params[position_of_ki]
        lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([ki >= 0.001 and ki < 1000]))

    return all(lie_in_correct_interval)

def count_neg_charged_residues(sequence):

    return sequence.count('D') + sequence.count('E')

def count_pos_charged_residues(sequence):

    return sequence.count('K') + sequence.count('R')

def aa_sequence_to_weight(sequence,monoisotopic=False):

    # https://proteomicsresource.washington.edu/protocols06/masses.php
    if not monoisotopic:

        amino_acid_weights = {
            "A": 71.0779, "R": 156.18568, "N": 114.10264, "D": 115.0874, "C": 103.1429,
            "Q": 128.12922, "E": 129.11398, "G": 57.05132, "H": 137.13928, "I": 113.15764,
            "L": 113.15764, "K": 128.17228, "M": 131.19606, "F": 147.17386, "P": 97.11518,
            "S": 87.0773, "T": 101.10388, "W": 186.2099, "Y": 163.17326, "V": 99.13106,
            "U": 150.3079,"O":237.29816
        }

    else:

        amino_acid_weights = {
            "A": 71.037113805, "R": 156.101111050, "N": 114.042927470, "D": 115.026943065, "C": 103.009184505,
            "Q": 128.058577540, "E": 129.042593135, "G": 57.021463735, "H": 137.058911875, "I": 113.084064015,
            "L": 113.084064015, "K": 128.094963050, "M": 131.040484645	, "F": 147.068413945, "P": 97.052763875,
            "S": 87.032028435, "T": 101.047678505, "W": 186.079312980, "Y": 163.063328575, "V": 99.068413945,
            "U": 150.953633405,"O":237.147726925
        }

    weight = 0.0
    for aa in sequence:
        if aa in amino_acid_weights:
            weight += amino_acid_weights[aa]
        else:
            continue

    # Add the N-terminus (H) and C-terminus (OH) groups to calculate the neutral mass of the peptide/protein.
    if monoisotopic:
        oxygen_mass = 15.99491463
        hydrogen_mass = 1.007825
    else:
        oxygen_mass   = 15.9994
        hydrogen_mass = 1.00794

    weight += oxygen_mass + 2 * hydrogen_mass

    # Round to the first digit
    weight = np.round(weight,1)

    return weight

def calculate_epsilon(sequence):
    """
    Given a protein sequence with natural aminoacids
    compute the molar extinction coefficient at 205 nm, 214 nm and 280 nm
    """

    # Count the number of times we have each aminoacid in the sequence string
    n_A = sequence.count('A')
    n_C = sequence.count('C')
    n_D = sequence.count('D')
    n_E = sequence.count('E')
    n_F = sequence.count('F')
    n_G = sequence.count('G')
    n_H = sequence.count('H')
    n_I = sequence.count('I')
    n_K = sequence.count('K')
    n_L = sequence.count('L')
    n_M = sequence.count('M')
    n_N = sequence.count('N')
    n_P = sequence.count('P')
    n_Q = sequence.count('Q')
    n_R = sequence.count('R')
    n_S = sequence.count('S')
    n_T = sequence.count('T')
    n_V = sequence.count('V')
    n_W = sequence.count('W')
    n_Y = sequence.count('Y')

    # Start with Method 214 nm - https://pubs.acs.org/doi/epdf/10.1021/jf070337l?ref=article_openPDF
    # Does not depend on the number of disulfide bonds

    # Check if we have a proline at the N-terminus
    proline_N = sequence[0] == 'P'

    if proline_N:
        n_P -= 1
        proline_N_contribution = 30
    else:
        proline_N_contribution = 0

    epsilon_214 = (
            proline_N_contribution + n_A * 32 + n_C * 225 +
            n_D * 58 + n_E * 78 + n_F * 5200 + n_G * 21 + n_H * 5125 +
            n_I * 45 + n_K * 41 + n_L * 45 + n_M * 980 + n_N * 136 +
            n_P * 2675 + n_Q * 142 + n_R * 102 + n_S * 34 + n_T * 41 +
            n_V * 43 + n_W * 29050 + n_Y * 5375 + (len(sequence)-1)*923
                   )

    molecular_weight = aa_sequence_to_weight(sequence)

    # Compute 0.1 absorbance at 280 nm
    abs01_214 = np.round(epsilon_214 / molecular_weight,2)

    # Count how many possible disulfide bonds we have - number of cysteines divided by 2
    disulfide_bonds = sequence.count('C') // 2

    epsilon_280 = []
    epsilon_205 = []

    possible_disulfide_bonds = [x for x in range(disulfide_bonds+1)]

    # Iterate over the number of disulfide bonds
    for n_DB in possible_disulfide_bonds:

        # Method 280 nm
        epsilon_280.append(1490 * n_Y + 5500 * n_W + n_DB * 125)

        # Method 205 nm - https://pubmed.ncbi.nlm.nih.gov/23526461/

        epsilon_205.append(
            n_W * 20400 + n_F * 8600 + n_Y * 6080 + n_H * 5200 + n_M * 1830 +
            n_R * 1350 + n_C * 690 + n_N * 400 + n_Q * 400 +
            (len(sequence)-1) * 2780 + 820 * n_DB
        )

    abs01_280 = [np.round(x / molecular_weight,2) for x in epsilon_280]
    abs01_205 = [np.round(x / molecular_weight,2) for x in epsilon_205]

    # Create a dataframe with four columns, epsilon 280, epsilon 214, epsilon 205 and number of disulfide bonds
    df1 = pd.DataFrame({
        'Epsilon_280nm [1/(M*cm)]': epsilon_280,
        'Epsilon_214nm [1/(M*cm)]': epsilon_214,
        'Epsilon_205nm [1/(M*cm)]': epsilon_205,
        '#Disulfide bonds': possible_disulfide_bonds})

    # Create a dataframe with four columns, epsilon 280, epsilon 214, epsilon 205 and number of disulfide bonds
    df2 = pd.DataFrame({
        'Abs(0.1%)_280nm (=1 g/l)': abs01_280,
        'Abs(0.1%)_214nm (=1 g/l)': abs01_214,
        'Abs(0.1%)_205nm (=1 g/l)': abs01_205,
        '#Disulfide bonds': possible_disulfide_bonds})

    return [df1, df2]