# import 
import os
from typing import Annotated, List, Dict, Tuple, Optional, Union, Any
## 3rd-party
from pydantic import BaseModel
from langgraph.prebuilt import create_react_agent
from langchain_core.messages import BaseMessage, HumanMessage, AIMessage, ToolMessage
from langchain_openai import ChatOpenAI
from langchain_core.tools import tool
from langchain_core.runnables.config import RunnableConfig

def create_planner_tool():
    model = ChatOpenAI(model="gpt-4o-mini", temperature=0.2)

    class Step(BaseModel):
        explanation: str
        output: str

    class Planning(BaseModel):
        steps: list[Step]

    @tool
    def invoke_planner(
        message: Annotated[str, "Message to the planner"],
    ) -> Annotated[str, "Response from the planner"]:
        """
        The planner will help you plan how to accomplish a task.
        Be specific about the task you need help with.
        """
        prompt = "\n".join([
            "You are an expert planner.",
            "Think step by step to help the user accomplish the task.",
            "The task: ",
            message
        ])
        response = model.with_structured_output(Planning, strict=True).invoke(prompt)
        steps = "\n".join([f" - {x.output}" for x in response.steps])
        return {
            "messages": [HumanMessage(content=steps, name="planner")]
        }
    return invoke_planner

def create_critic_tool():
    model = ChatOpenAI(model="gpt-4o-mini", temperature=0.2)

    class Reasons(BaseModel):
        explanation: str
        output: str

    class Critque(BaseModel):
        critique: str
        reasoning: list[Reasons]

    @tool
    def invoke_critic(
        message: Annotated[str, "Message to the critic"],
    ) -> Annotated[str, "Response from the critic"]:
        """
        The critic will help you critique a plan to accomplish a task.
        """
        prompt = "\n".join([
            "You are an expert critic.",
            "Think step by step to critique the plan.",
            message
        ])
        response = model.with_structured_output(Critque, strict=True).invoke(prompt)
        reasoning = "\n".join([f" - {x.output}" for x in response.reasoning])
        reasoning = f"{response.critique}\nMy reasoning:\n{reasoning}"
        return {
            "messages": [HumanMessage(content=reasoning, name="critic")]
        }
    
    return invoke_critic