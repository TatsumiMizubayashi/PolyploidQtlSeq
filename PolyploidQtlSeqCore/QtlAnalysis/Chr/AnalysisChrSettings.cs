namespace PolyploidQtlSeqCore.QtlAnalysis.Chr
{
    /// <summary>
    /// 解析対象染色体設定
    /// </summary>
    internal class AnalysisChrSettings
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
        /// 解析対象染色体設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public AnalysisChrSettings(IAnalysisChrSettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            ChrSizeThreshold = new ChrSizeThreshold(settingValue.ChrSizeThreshold, parameterDictionary, userOptionDictionary);
            AnalysisChrNames = new AnalysisChrNames(settingValue.AnalysisChrNames, parameterDictionary, userOptionDictionary);
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
