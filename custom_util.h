void cl_profile_print(cl_event *, int event_num, double *sum_list, int *times);
void get_event_time(struct timeval *, struct timeval *, double *, int *);
void print_kernel_execution_time(double *, int *, int);
void print_event_execution_time(double *, int *);
