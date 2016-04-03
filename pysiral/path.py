import os
import glob


def validate_directory(folder):
    """
    Check if folder str is a existing folder, creates the folder if ``folder``
    does not exist and is of valid notation
    Returns status flag (True: valid and existing, False: invalid)
    """
    if os.path.exists(folder):
        # Existing folder -> Valid
        return True
    if os.path.isabs(folder):
        # Valid notation, try to create folder
        try:
            os.mkdir(folder)
        except:
            # Failed to create folder -> Invalid
            return False
        # Folder created -> Valid
        return True
    # not existing, not valid notation -> Invalid
    return False


def file_basename(filename, fullpath=False):
    """
    Returns the filename without file extension of a give filename (or path)
    """
    strarr = os.path.split(filename)
    file_name = strarr[-1]
    basename = file_name.split(".")[0]
    if fullpath:
        basename = os.path.join(strarr[0], basename)
    # XXX: Sketchy, needs better solution (with access to os documentation)
    return basename


def get_filenames(directory, file_extension):
    """"
    Parses a directory for certain file extensions and returns a sorted list

    Arguments:
        directory: str
            Directory to be searched
        file_extension: str:
            Identifier for file type (e.g. '*.dat')
    """
    filenames = glob.glob(os.path.join(directory, file_extension))
    return sorted(filenames)


def filename_from_path(path):
    return os.path.split(path)[1]


def folder_from_filename(filename):
    return os.path.split(filename)[0]


def filename_from_path(path):
    return os.path.split(path)[1]
