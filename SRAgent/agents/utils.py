from typing import List, Dict, Any
from langchain_openai import ChatOpenAI
from langchain_core.prompts import PromptTemplate


def create_step_summary_chain(model: str="gpt-4o-mini", max_tokens: int=45):
    """
    Create a chain of tools to summarize each step in a workflow.
    Args:
        model: The OpenAI model to use for the language model.
        max_tokens: The maximum number of tokens to use for the summary.
    Returns:
        A chain of tools to summarize each step in a workflow.
    """
    # Create the prompt template
    template = "\n".join([
        "Concisely summarize the provided step in the langgraph workflow.",
        f"The summary must be {max_tokens} tokens or less.",
        "Do not use introductory words such as \"The workflow step involves\"",
        "Write your output as plain text instead of markdown.",
        "#-- The workflow step --#",
        "{step}"
    ])
    prompt = PromptTemplate(input_variables=["step"], template=template)

    # Initialize the language model
    llm = ChatOpenAI(model_name=model, temperature=0, max_tokens=max_tokens)

    # Return the LLM chain
    return prompt | llm
