"""Run coveralls only on travis."""
import os
import subprocess

import click


def echo_call(cmd):
    """Echo command."""
    click.echo(f"calling: {' '.join(cmd)}", err=True)


@click.command()
def main():
    """Run coveralls only on travis."""
    if os.getenv('TRAVIS'):
        cmd = ['coveralls']
        echo_call(cmd)
        subprocess.call(cmd)


if __name__ == '__main__':
    main()
