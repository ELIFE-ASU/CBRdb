import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor


def mp_calc(func, arg, n=mp.cpu_count()):
    """
    Executes a function in parallel using a process pool.

    Parameters:
    func (callable): The function to execute.
    arg (iterable): An iterable of arguments to pass to the function.
    n (int, optional): The number of worker processes to use. Default is the number of CPU cores.

    Returns:
    list: A list of results from the function executions.
    """
    with mp.Pool(n) as pool:
        results = pool.map(func, arg)
    return results


def mp_calc_star(func, args, n=mp.cpu_count()):
    """
    Executes a function in parallel using a process pool with multiple arguments.

    Parameters:
    func (callable): The function to execute.
    args (iterable): An iterable of argument tuples to pass to the function.
    n (int, optional): The number of worker processes to use. Default is the number of CPU cores.

    Returns:
    list: A list of results from the function executions.
    """
    with mp.Pool(n) as pool:
        results = pool.starmap(func, args)
    return results


def tp_calc(func, arg, n=mp.cpu_count()):
    """
    Executes a function in parallel using a thread pool.

    Works best for I/O-bound tasks.

    Parameters:
    func (callable): The function to execute.
    arg (iterable): An iterable of arguments to pass to the function.
    n (int, optional): The number of worker threads to use. Default is the number of CPU cores.

    Returns:
    list: A list of results from the function executions.
    """
    with ThreadPoolExecutor(max_workers=n) as executor:
        results = executor.map(func, arg)
    return results
