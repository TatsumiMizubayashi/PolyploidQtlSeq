﻿using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 塩基クオリティ
    /// </summary>
    public class BaseQuality
    {
        /// <summary>
        /// 塩基クオリティの最小値
        /// </summary>
        [Obsolete("削除予定")]
        public const int MINIMUM = 1;

        /// <summary>
        /// 塩基クオリティの最大値
        /// </summary>
        [Obsolete("削除予定")]
        public const int MAXIMUM = 50;

        /// <summary>
        /// 塩基クオリティの規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 15;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "q";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "quality";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Threshold for base quality. qualified_quality_phred in fastp.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -q option must be an integer greater than or equal to 1 and less than or equal to 50.";

        /// <summary>
        /// 塩基クオリティを作成する。
        /// </summary>
        /// <param name="quality">塩基クオリティ</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public BaseQuality(int quality, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, quality, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// 塩基クオリティを取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        internal string ToFastpArg()
        {
            return $"-q {Value}";
        }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }
    }
}
