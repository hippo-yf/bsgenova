import itertools
import time
from threading import Thread, Event


# def spin(msg: str, done: Event) -> None:
#     for char in itertools.cycle(r'\|/-'):
#         status = f'\r{char} {msg}'
#         print(status, end='', flush=True)
#         if done.wait(.1):
#             break
#     blanks = ' ' * len(status)
#     print(f'\r{blanks}\r', end='')
# def slow() -> int:
#     time.sleep(8)
#     return 42
# def supervisor() -> int:
#     done = Event()
#     spinner = Thread(target=spin, args=('thinking!', done))
#     print(f'spinner object: {spinner}')
#     spinner.start()
#     result = slow()
#     done.set()
#     spinner.join()
#     return result
# def main() -> None:
#     result = supervisor()
#     print(f'Answer: {result}')

import subprocess

if __name__ == '__main__':
        # subprocess.run(["dir", ">a"], shell=True)
        subprocess.run(["python", "./src/bs-snv-caller8.py", "-i", "./data/atcg.example", "-o", "data/sample", "-P", "4"], shell=True)
