# update the $PATH environment variable
Exec {
  path => [
		"/usr/local/sbin",
		"/usr/local/bin",
		"/usr/sbin",
		"/usr/bin",
		"/sbin",
		"/bin",
	]
}

# keep package information up to date
exec {
	"apt_update":
	command => "/usr/bin/apt-get update"
}

# install packages.
package {
	"wget":            ensure => installed, require => Exec ["apt_update"];
	"bzip2":           ensure => installed, require => Exec ["apt_update"];
	"tar":             ensure => installed, require => Exec ["apt_update"];
	"git":             ensure => installed, require => Exec ["apt_update"];
	"build-essential": ensure => installed, require => Exec ["apt_update"];
	"gzip":            ensure => installed, require => Exec ["apt_update"];
	"zlib1g-dev":      ensure => installed, require => Exec ["apt_update"];
	"ncurses-dev":     ensure => installed, require => Exec ["apt_update"];
}

# command line tasks
exec {

  # install python2.7
  'install_python':
    command   => 'sudo apt install python2.7',
    cwd       => '/usr/bin',
    creates   => '/usr/bin/python2.7';

  # install pip
  'install_pip':
    command   => 'sudo apt install python-pip',
    cwd       => '/usr/bin',
    creates   => '/usr/bin/pip'
    require   => Exec[ 'install_python' ];

  # install required python packages
  'install_packages':
    command   => 'pip install --user argparse numpy biopython datetime gzip pysam termcolor',
    cwd       => '/usr/bin',
    require   => Exec[ 'install_pip' ];

  # export PYTHONPATH
  'export_pythonpath':
    command   => 'echo "export PYTHONPATH=/home/ubuntu/bin:$PYTHONPATH" >> .bashrc',
    cwd       => '/home/ubuntu',
    require   => Exec[ 'install_packages' ];


  # clone nextPARS repository
#  'clone_nextPARS':
#    command   => 'git clone https://github.com/Gabaldonlab/nextPARS.git',
#    cwd       => '/home/vagrant',
#    creates   => '/home/vagrant/nextPARS',
#    require   => Package[ 'git' ];
#  'chown_nextPARS':
#    command   => 'chown -R vagrant:vagrant /home/vagrant/nextPARS',
#    cwd       => '/home/vagrant',
#    require   => Exec[ 'clone_nextPARS' ];
}
