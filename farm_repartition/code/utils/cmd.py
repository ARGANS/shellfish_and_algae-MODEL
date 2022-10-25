import subprocess

def run(cmd:str) -> None:
    process =  subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    print('C: ', cmd)
    while True:
        line = process.stdout.readline()
        if not line:
            break
        output = line.decode().rstrip()
        print(output)