namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 解析対象染色体名
    /// </summary>
    internal class AnalysisChrNames
    {
        private static readonly char[] _delimiter = new[] { ',' };

        /// <summary>
        /// 解析対象染色体名を作成する。
        /// </summary>
        /// <param name="value">カンマ区切り染色体名</param>
        public AnalysisChrNames(string value)
        {
            Value = value;
            Names = string.IsNullOrEmpty(Value)
                ? Array.Empty<string>()
                : Value.Split(_delimiter);
            HasNames = Names.Any();
        }

        /// <summary>
        /// オプション指定された値を取得する。
        /// </summary>
        internal string Value { get; }

        /// <summary>
        /// 解析対象染色体名を取得する。
        /// </summary>
        internal string[] Names { get; }

        /// <summary>
        /// 解析対象染色体名を持っているかどうかを取得する。
        /// </summary>
        internal bool HasNames { get; }

    }
}
