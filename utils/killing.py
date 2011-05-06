#Some older functionality for killing all child processes of pam. This is linux specific

    def kill_single_process(self,pid):
        """
        kill (or try to kill) single process
        """
        from common import get_pid_label
      
        while True:
            try:
                while True:
                    print_std_out_err("Going to kill (SIGTERM) dirac.x related process with PID="+str(pid))
                    print_std_out_err(get_pid_label(pid))
                    time.sleep(3)
                    os.kill(pid,signal.SIGTERM) # or signal.SIGKILL
                    # waint and verify
                    time.sleep(4)
                    killedpid, stat = os.waitpid(pid, os.WNOHANG)
                    print_std_out_err("killedpid,stat:"+killedpid+" "+stat)
                    if killedpid == 0:
                        print_std_out_err("dirac.x process(es),PID="+str(pid)+\
                        ", not killed, going to repeat shoot after few seconds")
                        time.sleep(3)
                    else:
                        print_std_out_err(" Yes ! This process was killed !")
                    break # leave the while cycle...
            except:
                break

    def kill_dirac_processes(self):
        """
        kill all dirac processes...
        """
        from common import is_windows
        if is_windows() : return 		
	
        self.get_dirac_pids()
        was_killed=True
      
        i=0
        while was_killed==True or i<=2:
            was_killed=False

            if self.parallel_run:
                if self.dirac_in_pid(self.dirac_p.pid):
                    self.kill_single_process(self.dirac_p.pid)
                    was_killed=True
                    i=0

            for pid in self.all_dirac_pids:
                if self.dirac_in_pid(pid):
                    self.kill_single_process(pid)
                    was_killed=True
                    i=0
            # get fresh load of dirac processes
            self.get_dirac_pids()
            i=i+1
            print 'i=',i

        self.print_all_processes()
        self.print_timeout_message()
        return signal.SIGTERM

    def get_dirac_pids(self):
        """
        get all dirac PIDs
        """
        # first, check if the current PID has not changed
        if ( self.current_pid != os.getpid() ):
            for stream in [sys.stdout,sys.stderr]:
                print  >> stream,'*** change in current pam script process PID ! '
                print  >> stream,'previous assigned PID:',self.current_pid,' current PID:',os.getpid()
                print  >> stream,'PID restored to current PID',os.getpid()
                stream.flush()
            # restore current PID
            self.current_pid=os.getpid()

        #  get the list of children's processes; use two different methods
        from common import get_all_children
        self.dirac_pids=get_all_children(self.current_pid)

        # miro: keep it ...
        self.get_pidof_pids()

        self.dirac_pids.append(self.dirac_p.pid)
        self.all_dirac_pids=set.union( set(self.dirac_pidof_pids),set(self.dirac_pids) )

        # print
        if not self.silent:
            for stream in [sys.stdout,sys.stderr]:
                print >> stream,\
                '         dirac.x subprocess PID:',self.dirac_p.pid
                print >> stream,\
                '  All dirac.x python found PIDs:',self.dirac_pids
                print  >> stream,\
                '   All dirac.x pidof found PIDs:',self.dirac_pidof_pids
                print  >> stream,\
                '       All unified dirac.x PIDs:',self.all_dirac_pids
                stream.flush()

    def print_all_processes(self):
        """
        print all remaining processes using system command "ps"
        TODO: adapt for Windows
        """
        sys.stdout.flush()
        sys.stderr.flush()
        # command - print all processes
        command="ps -fu "+self.user_name
        process_string=subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,\
                       stdout=subprocess.PIPE).communicate()[0]
        print_std_out_err("\n current processes left: \n--------------------------- \n"+process_string)

    def dirac_in_pid(self,pid):
        """
        checks if pid process contains dirac.x/launcher ,i.e. is dirac related (returns True), or not (returns False)
        TODO: adapt for Windows
        """
        command="ps -p "+str(pid)
        process_string=subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,\
                       stdout=subprocess.PIPE).communicate()[0]
        # check on presense of substrings
        if 'dirac.x' in process_string:
            return True
        from common import extract_path_file
        if self.parallel_run:
            #TODO: fix it here.... (FIX what??)
            z=extract_path_file()
            z.get_path_file(self.mpirun)
            if "mpirun" in process_string:
                return True
            if "poe" in process_string:
                return True
        return False 

    def get_pidof_pids(self):
        """
        get all dirac subproces PIDs using Linux's pidof command
        store them in the "self.dirac_pidof_pids" list
        TODO: adapt for Windows
        """
        command="pidof dirac.x"
        z = subprocess.Popen(command,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
        dirac_pidof_str=list( (z.communicate()[0]).split()  )
        # get PID of the parallel launcher...
        if self.parallel_run:
            z1=subprocess.Popen("pidof "+self.mpirun,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            dirac_pidof_parallel=(z1.communicate()[0]).split()
            dirac_pidof_str.extend(dirac_pidof_parallel)
        # turn all processes into int
        self.dirac_pidof_pids = [int(item) for item in dirac_pidof_str]
