namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Bcftools 変異検出オプション
    /// </summary>
    internal class BcfToolsVariantCallSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [MinmumBaseQuality.SHORT_NAME] = MinmumBaseQuality.LONG_NAME,
            [MinmumBaseQuality.LONG_NAME] = MinmumBaseQuality.LONG_NAME,

            [MinmumMappingQuality.SHORT_NAME] = MinmumMappingQuality.LONG_NAME,
            [MinmumMappingQuality.LONG_NAME] = MinmumMappingQuality.LONG_NAME,

            [AdjustMappingQuality.SHORT_NAME] = AdjustMappingQuality.LONG_NAME,
            [AdjustMappingQuality.LONG_NAME] = AdjustMappingQuality.LONG_NAME
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
        /// Bcftools変異検出オプションを作成する。
        /// </summary>
        /// <param name="optionValues">bcftools変異検出オプション値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public BcfToolsVariantCallSettings(IBcftoolsVariantCallSettingValue optionValues, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            MinmumBaseQuality = new MinmumBaseQuality(optionValues.MinBq, parameterDictionary, userOptionDictionary);
            MinmumMappingQuality = new MinmumMappingQuality(optionValues.MinMq, parameterDictionary, userOptionDictionary);
            AdjustMappingQuality = new AdjustMappingQuality(optionValues.AdjustMq, parameterDictionary, userOptionDictionary);
        }


        /// <summary>
        /// Base Quality最低値
        /// </summary>
        public MinmumBaseQuality MinmumBaseQuality { get; }

        /// <summary>
        /// Mapping Quality最低値
        /// </summary>
        public MinmumMappingQuality MinmumMappingQuality { get; }

        /// <summary>
        /// Adjust Mapping Quality
        /// </summary>
        public AdjustMappingQuality AdjustMappingQuality { get; }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                MinmumBaseQuality.ToParameterFileLine(),
                MinmumMappingQuality.ToParameterFileLine(),
                AdjustMappingQuality.ToParameterFileLine()
            };
        }
    }
}
