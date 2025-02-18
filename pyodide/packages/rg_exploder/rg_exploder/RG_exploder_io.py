from io import StringIO
import pyodide

JS=None

def is_local_file(filename):
    # TODO mravn -- check browser storage
    return True

def open_local_read(filename):
    # TODO mravn -- read from browser storage
    return pyodide.open_url(filename)

def open_local_write(filename):
    # TODO mravn -- store in browser storage
    return DevNullFile(filename)

def is_file(filename):
    return True

def open_read(filename):
    print('Reading hosted file ' + filename)
    return pyodide.open_url(filename)

def open_write(filename, buffer_size = 4194304):
    return BrowserFile(filename, buffer_size)

class DevNullFile:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def write(self, line):
        pass

    def close(self):
        pass

class BrowserFile:
    def __init__(self, name, buffer_size):
        self.name = name[name.rfind('/')+1:]
        self.buffer = []
        self.buffer_size = buffer_size
        self.bytes_in_buffer = 0
        self.bytes_in_total = 0
        self.short_lines_buffer = [] # short lines optimization
        self.bytes_in_short_lines_buffer = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def flush(self):
        if 0 < len(self.buffer):
            JS.postLines(self.name, self.buffer)
        self.buffer = []
        self.bytes_in_buffer = 0

    def close(self):
        self.flush_short_lines()
        self.flush()
        JS.postEndOfFile(self.name)
        print(str(self.bytes_in_total) + " bytes written to " + self.name)

    def write(self, line):
        bytes = len(line)
        self.bytes_in_total += bytes
        if self.buffer_size == 0:
            JS.postLine(self.name, line)
        elif bytes + self.bytes_in_short_lines_buffer < 4096:
            self.short_lines_buffer.append(line)
            self.bytes_in_short_lines_buffer += bytes
        else:
            self.flush_short_lines()
            self.buffer.append(line)
            self.bytes_in_buffer += bytes
            if (self.buffer_size < self.bytes_in_buffer):
                self.flush()

    def flush_short_lines(self):
        if 0 < self.bytes_in_short_lines_buffer:
            self.buffer.append("".join(self.short_lines_buffer))
            self.bytes_in_buffer += self.bytes_in_short_lines_buffer
            self.short_lines_buffer.clear()
            self.bytes_in_short_lines_buffer = 0
