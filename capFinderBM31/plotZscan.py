import capFinderBM31 as cf
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help = 'zscan file for the program to read', type = str)
    parser.add_argument('-s', '--nostdevs', help = 'number of standard deviations used to find peaks (default 1)', type = float, default = 1)
    parser.add_argument('-c','--capsize', type = float, default = 1, help = 'capillary size (will avoid finding 2 peaks within 1 capillary size) (default 1)')

    args = parser.parse_args()

    file = args.file
    noStdevs = args.nostdevs
    capsize = args.capsize

    x = cf.plotResults(file, noStdevs=noStdevs, capsize = capsize)
    print(x.tolist())

if __name__ ==  '__main__':
    main()
