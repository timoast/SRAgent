# import
## batteries
import os
import sys
import operator
from enum import Enum
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any, Sequence, TypedDict
## 3rd party
from dotenv import load_dotenv
from pydantic import BaseModel
from langchain_openai import ChatOpenAI
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage
from langgraph.graph import StateGraph, START, END
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder

# functions
class GraphState(TypedDict):
    """
    Shared state of nodes in the graph
    """
    messages: Annotated[Sequence[BaseMessage], operator.add]
    route: Annotated[str, "route to next node"]
    rounds: Annotated[int, operator.add, "number of rounds"]

## planner function
class Step(BaseModel):
    explanation: str
    output: str

class Planning(BaseModel):
    steps: list[Step]

def create_planner_node():
    model = ChatOpenAI(model="gpt-4o", temperature=0.2)

    def invoke_planner(
        state: GraphState
    ) -> Annotated[dict, "Response from the critic"]:
        """
        The planner tool takes in the current state of the conversation and generates a plan to achieve the goal.
        """
        # create prompt
        prompt = ChatPromptTemplate.from_messages([
            # First add any static system message if needed
            ("system", "You are an expert planner. Think step by step to critique the plan. Concisely explain your reasoning for each step in one or two sentences."),
            # Include all previous messages from the state
            MessagesPlaceholder(variable_name="history"),
            # Add the final question/instruction
            ("human", "Based on the messages above, create a new or revised plan to achieve the goal. The plan should be high-level and concise."),
        ])
        formatted_prompt = prompt.format_messages(
            history=state["messages"][-3:]
        )
        # call the model
        response = model.with_structured_output(Planning, strict=True).invoke(formatted_prompt)
        # format the response
        steps = "\n".join([f" - {x.output}" for x in response.steps])
        return {
            "messages": [AIMessage(content=steps, name="planner")]
        }

    return invoke_planner

## critic function
class Reasons(BaseModel):
    explanation: str
    output: str

class Critque(BaseModel):
    critique: str
    reasoning: list[Reasons]

def create_critic_node():
    model = ChatOpenAI(model="gpt-4o-mini", temperature=0.1)

    def invoke_critic(
        state: GraphState
    ) -> Annotated[dict, "Response from the critic"]:
        """
        The critic will help you critique a plan to accomplish a task.
        """
        # create prompt
        prompt = ChatPromptTemplate.from_messages([
            # First add any static system message if needed
            ("system", "You are an expert critic. Think step by step to critique the plan. Concisely explain your reasoning in one or two sentences."),
            # Include all previous messages from the state
            MessagesPlaceholder(variable_name="history"),
            # Add the final question/instruction
            ("human", "Based on the messages above, critique the current plan."),
        ])
        formatted_prompt = prompt.format_messages(
            history=state["messages"][-3:]
        )
        # call the model
        response = model.with_structured_output(Critque, strict=True).invoke(formatted_prompt)
        # format the response
        reasoning = "\n".join([f" - {x.output}" for x in response.reasoning])
        reasoning = f"{response.critique}\nMy reasoning:\n{reasoning}"
        return {
            "messages": [AIMessage(content=reasoning, name="critic")]
        }
    
    return invoke_critic

class Choices(Enum):
    CONTINUE = "Continue"
    STOP = "Stop"

class Choice(BaseModel):
    Choice: Choices

def create_router_node():
    model = ChatOpenAI(model="gpt-4o-mini", temperature=0)

    def invoke_router(
        state: GraphState
    ) -> Annotated[dict, "Response from the router"]:
        """
        Route the conversation to the appropriate tool based on the current state of the conversation.
        """
        # create prompt
        prompt = ChatPromptTemplate.from_messages([
            # First add any static system message if needed
            ("system", "You determine whether to continue or stop the conversation."),
            # Include all previous messages from the state
            ("system", "Here is the task:"),
            MessagesPlaceholder(variable_name="task"),
            ("system", "Here are the last few messages:"),
            MessagesPlaceholder(variable_name="history"),
            # Add the final question/instruction
            ("human", "Based on the messages above, select STOP if the task is complete or CONTINUE for another round of planning and critique."),
        ])
        # get the last messages, but it cannot include the first message
        last_n_messages = -6 if len(state["messages"]) > 6 else -(len(state["messages"]) - 1)
        formatted_prompt = prompt.format_messages(
            task=[state["messages"][0].content],
            history=state["messages"][last_n_messages:]
        )
        # call the model
        response = model.with_structured_output(Choice, strict=True).invoke(formatted_prompt)
        return {"route": response.Choice.value, "rounds": 1}

    return invoke_router

def route_interpret(state: GraphState) -> str:
    """
    Determine the route based on the current state of the conversation.
    """
    if state["rounds"] >= 2:
        return END
    return "planner_node" if state["route"] == "Continue" else END

def create_planner_graph():
    """
    Create a graph that combines the planner, critic, and router nodes.
    """
    workflow = StateGraph(GraphState)

    # nodes
    workflow.add_node("planner_node", create_planner_node())
    workflow.add_node("critic_node", create_critic_node())
    workflow.add_node("router_node", create_router_node())

    # edges
    workflow.add_edge(START, "planner_node")
    workflow.add_edge("planner_node", "critic_node")
    workflow.add_edge("critic_node", "router_node")
    workflow.add_conditional_edges("router_node", route_interpret)

    # compile the graph
    graph = workflow.compile()
    return graph

# main
if __name__ == "__main__":
    # setup
    from dotenv import load_dotenv
    load_dotenv()

    # planner
    state = {"messages" : [
        HumanMessage(content="How do I convert GSE121737 to an SRA accession?", name="human")
    ]}
    # invoke_planner(state)

    # critic
    state = {"messages" : [
        HumanMessage(content="How do I convert GSE121737 to an SRA accession?", name="human"),
        AIMessage(content="Run Bio.Entrez.esearch ", name="planner"),
    ]}
    # invoke_critic(state)

    # test
    state = {"messages" : [
        HumanMessage(content="How do I convert GSE121737 to an SRA accession?", name="human"),
        AIMessage(content="Run Bio.Entrez.esearch ", name="planner"),
        AIMessage(content="The plan is good! You can finish", name="critic"),
    ]}
    #invoke_router(state)

    graph = create_planner_graph()

    # Call the graph
    input = {"messages" : [
        HumanMessage(content="How do I convert GSE121737 to an SRA accession?", name="human")
    ]}
    for step in graph.stream(input, config={"max_concurrency" : 2, "recursion_limit": 30}):
        print(step)