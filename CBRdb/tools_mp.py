import fcntl
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable, Any, Tuple


def mp_calc(func: Callable[[Any], Any], arg: Iterable[Any], n: int = mp.cpu_count()) -> list[Any]:
    """
    Executes a function in parallel using a process pool.

    Parameters
    ----------
    func : Callable[[Any], Any]
        The function to execute.
    arg : Iterable[Any]
        An iterable of arguments to pass to the function.
    n : int, optional
        The number of worker processes to use. Default is the number of CPU cores.

    Returns
    -------
    list[Any]
        A list of results from the function executions.
    """
    with mp.Pool(n) as pool:
        results = pool.map(func, arg)
    return results


def mp_calc_star(func: Callable[..., Any], args: Iterable[Tuple[Any, ...]], n: int = mp.cpu_count()) -> list[Any]:
    """
    Executes a function in parallel using a process pool with multiple arguments.

    Parameters
    ----------
    func : Callable[..., Any]
        The function to execute.
    args : Iterable[Tuple[Any, ...]]
        An iterable of argument tuples to pass to the function.
    n : int, optional
        The number of worker processes to use. Default is the number of CPU cores.

    Returns
    -------
    list[Any]
        A list of results from the function executions.
    """
    with mp.Pool(n) as pool:
        results = pool.starmap(func, args)
    return results


def tp_calc(func: Callable[[Any], Any], arg: Iterable[Any], n: int = mp.cpu_count()) -> list[Any]:
    """
    Executes a function in parallel using a thread pool.

    Works best for I/O-bound tasks.

    Parameters
    ----------
    func : Callable[[Any], Any]
        The function to execute.
    arg : Iterable[Any]
        An iterable of arguments to pass to the function.
    n : int, optional
        The number of worker threads to use. Default is the number of CPU cores.

    Returns
    -------
    list[Any]
        A list of results from the function executions.
    """
    with ThreadPoolExecutor(max_workers=n) as executor:
        results = list(executor.map(func, arg))
    return results


def mp_calc_chunked(
        func: Callable[[Any], Any],
        arg: Iterable[Any],
        n: int | None = None,
        chunksize: int | None = None,
) -> list[Any]:
    """
    Executes a function in parallel using a process pool, with optional chunking.

    Parameters
    ----------
    func : Callable[[Any], Any]
        The function to execute (on a single element).
    arg : Iterable[Any]
        An iterable of arguments to pass to the function.
    n : int, optional
        Number of worker processes (default: mp.cpu_count()).
    chunksize : int, optional
        How many items each worker gets per batch. If None, multiprocessing
        chooses a default based on len(arg).

    Returns
    -------
    list[Any]
        A list of results from the function executions.
    """
    if n is None:
        n = mp.cpu_count()

    with mp.Pool(n) as pool:
        results = pool.map(func, arg, chunksize=chunksize or 1)
    return results


def write_to_shared_file(message: str, shared_file: str) -> None:
    """
    Write a message to a shared file with an exclusive lock.

    Args:
        message (str): The message to write to the file.
        shared_file (str): The path to the shared file.

    Returns:
        None
    """
    with open(shared_file, 'a') as f:
        # Acquire an exclusive lock before writing
        fcntl.flock(f, fcntl.LOCK_EX)
        # Write the message to the file
        f.write(message)
        # Release the lock after writing
        fcntl.flock(f, fcntl.LOCK_UN)
    return None
