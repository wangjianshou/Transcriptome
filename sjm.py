class Topo(object):
  '''
  Used to build sjm job.
  '''
  def __init__(self, job):
    self.job = job
    self.joblist = []

  def initialize(self):
    self.joblist=[]

  def __call__(self, *args, **kwargs):
    self.joblist.append(self.job(*args, **kwargs))

  def __iter__(self):
    self.idx = 0
    self.order = []
    return self

  def __next__(self):
    if self.idx == len(self.joblist):
      raise StopIteration
    else:
      order = self.joblist[self.idx].order.strip()
      order != '' and self.order.append(order)
      result = self.joblist[self.idx].sjm
      self.idx += 1
      return result

@Topo
class Job(object):
  '''
  This class is used for building sjm job.
  '''
  __slot__ = ['name', 'ater', 'sched', 'sjm', 'number', '__dict__']

  def __init__(self, name, shell, after=None, vf=None, proc=None, local=False, project=None, queue=None):
    '''
    Type of all parameters should be string except after, unless you use default.
    Parameter of after can be string, tuple or list which element is string.
    '''
    self.name = name
    name = '\tname ' + name
    self.after = set()
    if isinstance(after, str):
      self.after.add(after)
    else:
      bool(after) and self.after.update(after)
    #self.number = len(self.after)
    project = project and '-P ' + project
    queue =  queue and '-q ' + queue
    try:
      sched = ' '.join(("\tsched_options", project, queue, "-cwd", "-l vf="+vf+",p="+proc))
    except TypeError:
      sched = ' '.join(("\tsched_options", "-V -cwd", "-l vf="+vf+",p="+proc))
    finally:
      sched = local and "\thost localhost" or sched
    cmd = "\tcmd sh " + shell
    self.sjm = '\n'.join(['job_begin', name, sched, cmd, 'job_end'])

  @property
  def order(self):
    eachorder = tuple(map(lambda x: 'order '+self.name+' after '+x, self.after))
    return '\n'.join(eachorder)

  '''def delete(self, jobname):
    self.number = jobname in self.after and self.number-1 or self.number
    if self.number == 0:
      return self
  '''
