# import
## batteries
import os
import sys
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any
## 3rd party
from dotenv import load_dotenv
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent
## package
from SRAgent.tools.planner import create_planner_tool, create_critic_tool

# functions
def create_planner_agent():
    model = ChatOpenAI(model="gpt-4o", temperature=0.1)

    state_mod = "\n".join([
        "You are a supervisor.",
        "Your job is to help the user plan and critique a plan to accomplish a task.",
        "Facilitate the conversation between the planner and the critic.",
        "When the planner is done, ask the critic to critique the plan.",
        "When the critic is done, ask the planner to revise the plan.",
        "Use up to 3 planner-critic rounds to accomplish the task.",
        "Return the final plan.",
    ])
    planner_agent = create_react_agent(
        model=model,
        tools=[create_planner_tool(), create_critic_tool()],
        state_modifier=state_mod
    )
    return planner_agent


# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()

    # create and run agent
    planner_agent = create_planner_agent()
    input = {"messages": [("user", "Convert GSE121737 to SRX accessions")]}
    print(planner_agent.invoke(input))