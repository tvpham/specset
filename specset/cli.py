import argparse


def hello_world():
    print("Hello, world")


def main():
    parser = argparse.ArgumentParser(prog="specset", description="Toolkit for manipulation of a set of mass spectra")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("hello", help="Print Hello, world")

    args = parser.parse_args()

    if args.command == "hello":
        hello_world()


if __name__ == "__main__":
    main()
