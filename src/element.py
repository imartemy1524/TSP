class Element:
    def __init__(self, to: int, count: int):
        self.to = to
        self.count = count

    def __str__(self):
        return f'{self.to} ({self.count})'
