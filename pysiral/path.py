import os
import glob


def validate_directory(dir):
    """
    Tests if given directory exists and creates directory if not

    Arguments:
        dir: str
            Path to the directory
    """
    if not os.path.isdir(dir):
        os.mkdir(dir)


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
