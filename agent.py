class Agent:
    # Constructor
    def __init__(self, x, y, state):
        self.x = x # x coordinate
        self.y = y # y coordinate
        self.state = state # Empty = 0 / S = 1 / E = 2 / I = 3 / R = 4
    
    """
    get_position and get_state are not as important here in Python
    because the attributes are already public but they would be useful
    if the code was written in c++ with private and public class
    """
    # Return agent's position
    def get_position(self):
        return self.x, self.y
    
    # Return agent's SEIR state
    def get_state(self):
        return self.state
    
    # Update the agent's position
    def set_position(self, xnew, ynew):
        self.x = xnew
        self.y = ynew

    # Update agent's SEIR state
    def set_state(self, new_state):
        self.state = new_state