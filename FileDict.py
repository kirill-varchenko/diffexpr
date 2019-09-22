import re


class FileDict(dict):
    """Associate dict with file and given regex pattern with 'key' and 'value' named groups"""
    def __init__(self, file, pattern):
        super().__init__()
        self.file = file
        self.pattern = pattern
        self.re_comp = re.compile(pattern)

    def preload(self, *keys):
        """Preload items with given keys from file to dict"""
        set_keys = {str(k) for k in keys}
        with open(self.file, 'r') as fi:
            for line in fi:
                m = self.re_comp.search(line)
                key, value = m.group('key'), m.group('value')
                if key in set_keys:
                    set_keys.remove(key)
                    self[key] = value
                if not set_keys:
                    break

    def __getitem__(self, item):
        if item not in self.keys():
            self.preload(item)
        return super().__getitem__(item)
