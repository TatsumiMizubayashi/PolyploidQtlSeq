using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Adjust Mapping Quality
    /// </summary>
    public class AdjustMappingQuality
    {
        /// <summary>
        /// adjust MQの最小値
        /// </summary>
        private const int MINIMUM = 0;

        /// <summary>
        /// adjust MQの最大値
        /// </summary>
        private const int MAXIMUM = 1000;

        /// <summary>
        /// adjust MQの規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 60;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "C";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "adjustMQ";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Value for adjust mapping quality at variant detection in bcftools mpileup. \r\nSpecify 0, to disable this function.\r\n";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -C option must be an integer greater than or equal to 0 and less than or equal to 100.";

        /// <summary>
        /// Adjust Mapping Qualityを作成する。
        /// </summary>
        /// <param name="adjustMq">adjustMq</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public AdjustMappingQuality(int adjustMq, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, adjustMq, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }


        /// <summary>
        /// Adjust Mapping Qualityを作成する。
        /// </summary>
        /// <param name="adjustMq">adjustMq</param>
        public AdjustMappingQuality(int adjustMq)
        {
            if (adjustMq < MINIMUM || adjustMq > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(adjustMq));

            Value = adjustMq;
        }

        /// <summary>
        /// Adjust Mapping Qualityを取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }

    }
}
