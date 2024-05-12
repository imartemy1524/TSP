SIZE = 120
data = '''

'''

data_list = [[int(j.strip()) for j in i.split(' ') if j] for i in data.split('\n')]
data_list_1d = [j for i in data_list for j in i]
print(len(data_list_1d)/120)
