

# based on https://www.w3resource.com/python-exercises/tkinter/python-tkinter-widgets-exercise-18.php
#def sort_descending(treeframe, column):
#    # first we need the data as a list of tupes row by row
#    data = []
#    values = treeframe.get_children('')
#    for node in values:
#        data.append(treeframe.set(node, column), node)        # treeframe.set(row, column) gets a value at that position
#    # we sort the data 
#    data.sort(reverse=True)
#    # chaos: let's mix it up by moving the items around
#    for index, (value, node) in enumerate(data):
#        treeframe.move(node, '', index)
#    treeframe.heading(column, command=lambda: sort_descending(treeframe, column))
# the following method did not work and I ran out of time to fix it

def create_headings(table, headings):
    # headings provided as a list of strings
    headingsCount = len(headings)
    count = 0
    for heading in headings:
        table.heading(headings[count], text = headings[count]) #command=lambda: sort_descending(table, headings[count-1]))
        count +=1

def TheChosen(table):
    chosen = table.selection()
    ID = []
    for entry in chosen:
        ID.append(table.item(entry)["values"])
    print(ID)