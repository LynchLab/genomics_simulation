import resource

# the soft limit imposed by the current configuration
# the hard limit imposed by the operating system.
soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
print 'Soft limit is ', soft, "and the hard limit is ", hard

# For the following line to run, you need to execute the Python script as root.
resource.setrlimit(resource.RLIMIT_NOFILE, (3000, hard))
