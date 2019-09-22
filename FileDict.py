import re


class FileDict(dict):
    """Associate dict with file and given regex pattern with 'key' and 'value' named groups"""

    def __init__(self, file, pattern, key_converter=None, value_converter=None):
        super().__init__()
        self.file = file
        self.pattern = pattern
        self.re_comp = re.compile(pattern)
        self.key_converter = key_converter
        self.value_converter = value_converter

    def preload(self, *keys):
        """Preload items with given keys from file to dict"""
        set_keys = set(keys)
        with open(self.file, 'r') as fi:
            for line in fi:
                m = self.re_comp.search(line)
                key = m.group('key') if self.key_converter is None else self.key_converter(m.group('key'))
                value = m.group('value') if self.value_converter is None else self.value_converter(m.group('value'))
                if key in set_keys:
                    set_keys.remove(key)
                    self[key] = value
                if not set_keys:
                    break

    def __getitem__(self, item):
        if item not in self.keys():
            self.preload(item)
        return super().__getitem__(item)
