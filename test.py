class MyDictionary:
    def __init__(self):
        self.list = []
        for i in range(100):
            self.list.append([])

    def get(self, key):
        index = self.GenerateHash(key)
        memo = self.list[index]
        for k, v in memo:
            if k == key:
                return v
        raise KeyError(f"Key {key} not found")

    def add(self, key, value):
        index = self.GenerateHash(key)
        memo = self.list[index]
        for i, (k, v) in enumerate(memo):
            if k == key:
                memo[i] = (k, value)
                return
        memo.append((key, value))

    def remove(self, key):
        index = self.GenerateHash(key)
        memo = self.list[index]
        for i, (k, v) in enumerate(memo):
            if k == key:
                del memo[i]
                return
        raise KeyError(f"Key {key} not found")

    def GenerateHash(self, key):
        return hash(key) % len(self.list)


my_dict = MyDictionary()
my_dict.add("key", "value")
print(my_dict.get("key"))
my_dict.remove("key")
try:
    my_dict.get("key")
except KeyError as e:
    print(e)
