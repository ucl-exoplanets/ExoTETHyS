import os, click, shutil


def ls_database(grid='', starts=None, ends=None, contains=None):
    """
    By default this function lists the path and content of the ExoTETHyS database. It can be used to list all the files or a subset of files from any stellar models grid directory.
    
    :argument str grid: name of the stellar models grid (default is empty string)
    :argument str starts: to list the name of the files starting with the given string, for the selected grid (default is None)
    :argument str ends: to list the name of the files starting with the given string, for the selected grid (default is None)
    :argument list of str contains: to list the name of the files containing the given strings (default is None)
    ..note:: the arguments starts, ends and contains are ignored if the grid is not specified
    :return: by default, the database path and the names of the folders in it, or the requested grid path and files in it
    :rtype: str, list of str
    """
    database_path = os.path.join(os.path.expanduser('~'), '.{0}'.format("exotethys"))
    path = os.path.join(database_path, grid)
    files = os.listdir(path)
    if starts:
        files = [f for f in files if f.startswith(starts)]
    if ends:
        files = [f for f in files if f.endswith(ends)]
    if contains:
        for piece in contains:
            files = [f for f in files if piece in f]
    return path, files


def cp_database(dir_item, final_path, files=None):
    """
    This function can copy (part of) the ExoTETHyS database to another destination.
    
    :param str dir_item: it can be 'all' or the name of a stellar models grid
    :param str final_path: path of the existing directory where to copy the requested items
    :argument list of str files: names of the files from the selected grid to copy (default is None)
    :return:
    """
    if not os.path.exists(final_path):
        print('ERROR: the desired final path does not exists.')
    else:
        [database_path, grids] = ls_database()
        if dir_item == 'all':
            very_final_path = os.path.join(final_path, 'ExoTETHyS_database')
            shutil.copytree(database_path, very_final_path)
        elif dir_item in grids:
            grid_to_copy = os.path.join(database_path, dir_item)
            if not files:
                very_final_path = os.path.join(final_path, dir_item)
                shutil.copytree(grid_to_copy, very_final_path)
            else:
                for file in files:
                    file_to_copy = os.path.join(grid_to_copy, file)
                    very_final_path = os.path.join(final_path, file)
                    shutil.copy(file_to_copy, very_final_path)
        else:
            print('ERROR:', dir_item, 'unrecognised item to copy.')


def rm_database(grid=None, files=None):
    """
    This function can remove (part of) the ExoTETHyS database. By default, it removes the whole database. The action is irreversible. For safety, the user will be asked to confirm their intention.
    
    :argument str grid: the name of the stellar models grid to remove or containing the files to remove (default is None, to remove the whole database)
    :argument list of str files: names of the files from the selected grid to remove (default is None, to remove the whole grid)
    :return:
    """
    database_path = os.path.join(os.path.expanduser('~'), '.{0}'.format("exotethys"))
    #to remove specific files
    if files:
        if not grid:
            print('ERROR: grid must be specified to remove files')
            return
        else:
            path = os.path.join(database_path, grid)
            for file in files:
                file_to_delete = os.path.join(path, file)
                message = 'Are you sure that you want to delete the file ' + file_to_delete + '?'
                sure = click.confirm(message)
                if sure:
                    if os.path.exists(file_to_delete):
                        os.remove(file_to_delete)
                    else:
                        print('ERROR: The file', file_to_delete, 'does not exists.')
        #return
    #to remove a full grid of files
    elif grid:
        path = os.path.join(database_path, grid)
        message = 'Are you sure that you want to delete the directory ' + path + '?'
        sure = click.confirm(message)
        if sure:
            if os.path.exists(path):
                shutil.rmtree(path)
            else:
                print('ERROR: The directory', path, 'does not exists.')
    #to remove the whole database
    else:
        message = 'Are you sure that you want to delete the directory ' + database_path + '?'
        sure = click.confirm(message, abort=True)
        if sure:
            if os.path.exists(database_path):
                shutil.rmtree(database_path)
            else:
                print('ERROR: The directory', path, 'does not exists.')
    #return ??????
