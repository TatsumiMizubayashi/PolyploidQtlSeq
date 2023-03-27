namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 解析対象染色体オプション
    /// </summary>
    internal class AnalysisChrOption
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [ChrSizeThreshold.SHORT_NAME] = ChrSizeThreshold.LONG_NAME,
            [ChrSizeThreshold.LONG_NAME] = ChrSizeThreshold.LONG_NAME,

            [AnalysisChrNames.SHORT_NAME] = AnalysisChrNames.LONG_NAME,
            [AnalysisChrNames.LONG_NAME] = AnalysisChrNames.LONG_NAME
        };

        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }

        /// <summary>
        /// 解析対象染色体オプションを作成する。
        /// </summary>
        /// <param name="optionValues">オプション値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public AnalysisChrOption(IAnalysisChrOptions optionValues, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            ChrSizeThreshold = new ChrSizeThreshold(optionValues.ChrSizeThreshold, parameterDictionary, userOptionDictionary);
            AnalysisChrNames = new AnalysisChrNames(optionValues.AnalysisChrNames, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// 染色体サイズしきい値
        /// </summary>
        public ChrSizeThreshold ChrSizeThreshold { get; }

        /// <summary>
        /// 解析対象染色体名
        /// </summary>
        public AnalysisChrNames AnalysisChrNames { get; }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                ChrSizeThreshold.ToParameterFileLine(),
                AnalysisChrNames.ToParameterFileLine(),
            };
        }
    }
}
